"""The novelty check itself: look-up, prediction, comparison, verdict."""
from datetime import datetime, timezone

from rdkit import Chem

from . import __version__, solvents
from .chem import identify, parse_reaction, parse_smiles
from .match import OVERLAP, agreement, align, centres, overlaps
from .parsers import norm_solvent
from .predict import predict
from .sources import SOURCES

NUCLEI = ("13C", "1H")
AGREE = 0.5          # agreement score from which a record counts as consistent
MAX_LISTED = 10      # records listed per nucleus in the report
SOURCE_RANK = {"nmrshiftdb2": 0, "nmrexp": 1, "nmrbank": 2}

VERDICT_TEXT = {
    "reported": "Experimental NMR data for this exact structure are present "
                "in the indexed public sources.",
    "reported_unverified": "Records for this exact structure exist, but none "
                           "of their peak lists fits the predicted spectrum; "
                           "the structure-spectrum pairing may be wrong.",
    "reported_as_variant": "No record for this exact structure, but the same "
                           "connectivity (other stereo or isotope form) has "
                           "experimental data.",
    "listed_without_experiment": "The structure is listed, but only with "
                                 "calculated spectra.",
    "not_found": "No record for this structure in the indexed sources. This "
                 "is not proof that it is unpublished.",
    "invalid_input": "The structure could not be parsed.",
}

LIMITATIONS = [
    "Coverage is limited to the indexed sources. NMRexp is drawn from "
    "supporting information of 2010-2024 papers; compounds published "
    "elsewhere, in patents, or in subscription databases (SciFinder, Reaxys) "
    "are not seen.",
    "NMRexp was extracted automatically from PDFs. Structures can be paired "
    "with the wrong paragraph and peak lists can contain recognition errors.",
    "Structures are compared by InChIKey. Salts, solvates, tautomers and "
    "different charge states are treated as different compounds.",
    "Predicted shifts come from an environment look-up trained on "
    "nmrshiftdb2. Solvent, concentration, temperature and stereochemistry "
    "are ignored, and exchangeable protons are unreliable.",
    "Peak assignment is by nearest predicted shift. It is a proposal, not a "
    "2D-NMR-grade assignment.",
    "No user spectrum is involved: the tool checks whether data have been "
    "published, it does not verify a sample.",
]


class Log:
    """Collects every place where the run fell back to something weaker."""

    def __init__(self):
        self.entries = []

    def __call__(self, code, subject, detail):
        entry = {"code": code, "subject": subject, "detail": detail}
        if entry not in self.entries:
            self.entries.append(entry)


def _record_view(rec, match, score):
    return {
        "evidence": rec["origin"],
        "match": match,
        "source": rec["source"],
        "record_id": rec["record_id"],
        "doi": rec["doi"] or None,
        "nucleus": rec["nucleus"],
        "solvent": norm_solvent(rec["solvent"]) or None,
        "frequency": rec["frequency"] or None,
        "peaks": rec["peaks"],
        "agreement": score,
    }


def _labelled_smiles(mol):
    """SMILES with atom indices as map numbers, so assignments can be read."""
    tagged = Chem.Mol(mol)
    for atom in tagged.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    return Chem.MolToSmiles(tagged)


def _atom_label(mol, atoms):
    symbol = mol.GetAtomWithIdx(atoms[0]).GetSymbol()
    return symbol + ",".join(str(a + 1) for a in atoms)


def describe(index, smiles, role, log):
    """Everything known about one species."""
    out = {"role": role, "input": smiles}
    mol = parse_smiles(smiles)
    ident = identify(mol) if mol is not None else None
    if ident is None:
        log("invalid_structure", smiles, "RDKit could not parse this SMILES")
        out["status"] = "invalid_input"
        return out, None

    out.update(smiles=ident.smiles, inchikey=ident.inchikey,
               formula=ident.formula, atom_map=_labelled_smiles(mol))

    exact = index.compounds(inchikey=ident.inchikey)
    variants = [c for c in index.compounds(skeleton=ident.skeleton)
                if c["inchikey"] != ident.inchikey]
    exact_ids = [c["id"] for c in exact]
    variant_ids = [c["id"] for c in variants]

    predicted, spectra, best = {}, {}, {}
    for nucleus in NUCLEI:
        pred = predict(mol, nucleus, index.lookup_env)
        predicted[nucleus] = pred
        weak = [p for p in pred if p["radius"] < 2]
        if any(p["radius"] < 0 for p in pred):
            log("prediction_fallback", f"{ident.smiles} {nucleus}",
                "some atoms have no trained environment; fixed atom-type "
                "values used")
        elif weak:
            log("prediction_low_radius", f"{ident.smiles} {nucleus}",
                f"{len(weak)} of {len(pred)} signals predicted from an "
                "environment of radius < 2 bonds")
        shifts = [p["shift"] for p in pred]

        views = []
        for match, ids in (("exact", exact_ids), ("variant", variant_ids)):
            for rec in index.spectra(ids, nucleus):
                score = agreement(nucleus, shifts, centres(nucleus, rec["peaks"]))
                views.append(_record_view(rec, match, score))
        views.sort(key=lambda v: (v["match"] != "exact",
                                  v["evidence"] != "experimental",
                                  -v["agreement"]["score"],
                                  SOURCE_RANK.get(v["source"], 9)))
        spectra[nucleus] = views
        usable = [v for v in views
                  if v["match"] == "exact" and v["evidence"] == "experimental"]
        best[nucleus] = usable[0] if usable else None

    out["predicted"] = predicted
    out["records"] = spectra
    out["reference"] = best
    in_training = any(v["source"] == "nmrshiftdb2" and v["match"] == "exact"
                      for n in NUCLEI for v in spectra[n])
    out["prediction_independent"] = not in_training
    if in_training:
        log("prediction_not_independent", ident.smiles,
            "this compound is part of the predictor's training data, so "
            "agreement with nmrshiftdb2 records is partly circular")
    return out, mol


def _verdict(index, info):
    views = [v for n in NUCLEI for v in info["records"][n]]
    exact = [v for v in views if v["match"] == "exact"]
    exact_exp = [v for v in exact if v["evidence"] == "experimental"]
    variant_exp = [v for v in views
                   if v["match"] == "variant" and v["evidence"] == "experimental"]
    required = [n for n, s in SOURCES.items() if s.required]
    loaded = [n for n in required if n in index.meta["sources"]]
    coverage = len(loaded) / len(required)
    reasons = []

    def consistent(items):
        return [v for v in items if v["agreement"]["score"] >= AGREE]

    if exact_exp:
        docs = {v["doi"] or f'{v["source"]}:{v["record_id"].split("/")[0]}'
                for v in exact_exp}
        srcs = {v["source"] for v in exact_exp}
        good = consistent(exact_exp)
        conf = 0.70
        reasons.append(f"{len(exact_exp)} experimental record(s) with an "
                       "identical InChIKey (base 0.70)")
        if good:
            conf += 0.10
            reasons.append(f"{len(good)} record(s) fit the predicted spectrum "
                           "(+0.10)")
        if len(docs) >= 2:
            conf += 0.10
            reasons.append(f"{len(docs)} independent documents (+0.10)")
        if len(srcs) >= 2:
            conf += 0.05
            reasons.append("found in more than one database (+0.05)")
        status = "reported"
        if not good:
            conf -= 0.25
            status = "reported_unverified"
            reasons.append("no record fits the predicted spectrum (-0.25)")
        return status, round(min(conf, 0.95), 2), reasons

    if variant_exp:
        conf = 0.50
        reasons.append("experimental record(s) share the connectivity block of "
                       "the InChIKey only (base 0.50)")
        if consistent(variant_exp):
            conf += 0.10
            reasons.append("their peak lists fit the predicted spectrum (+0.10)")
        return "reported_as_variant", round(conf, 2), reasons

    if exact:
        return ("listed_without_experiment", 0.50,
                ["only calculated spectra are stored for this structure"])

    conf = round(0.5 * coverage, 2)
    reasons.append(f"absence can only be asserted for the indexed sources; "
                   f"{len(loaded)} of {len(required)} core sources loaded "
                   f"(0.50 x {coverage:.2f})")
    return "not_found", conf, reasons


def _assignment(info, mol, nucleus):
    pred = info["predicted"][nucleus]
    ref = info["reference"][nucleus]
    rows = []
    if not pred:
        return rows
    if ref is None:
        for p in pred:
            rows.append({"atoms": _atom_label(mol, p["atoms"]), "observed": None,
                         "predicted": p["shift"], "delta": None,
                         "radius": p["radius"], "support": p["support"],
                         "count": p["count"]})
        return rows
    observed = centres(nucleus, ref["peaks"])
    pairs, unexplained = align([p["shift"] for p in pred], observed)
    for i, j in pairs:
        p = pred[i]
        rows.append({"atoms": _atom_label(mol, p["atoms"]),
                     "observed": observed[j], "predicted": p["shift"],
                     "delta": round(observed[j] - p["shift"], 2),
                     "radius": p["radius"], "support": p["support"],
                     "count": p["count"]})
    for j in unexplained:
        rows.append({"atoms": None, "observed": observed[j], "predicted": None,
                     "delta": None, "radius": None, "support": None,
                     "count": None})
    rows.sort(key=lambda r: -(r["observed"] if r["observed"] is not None
                              else r["predicted"]))
    return rows


def _neighbours(index, info, limit=5):
    """Known isomers whose 13C spectrum is closest to the predicted one."""
    shifts = [p["shift"] for p in info["predicted"]["13C"]]
    if not shifts:
        return []
    isomers = {c["id"]: c for c in index.compounds(formula=info["formula"])
               if c["inchikey"] != info["inchikey"]}
    best = {}
    for rec in index.spectra(list(isomers), "13C"):
        if rec["origin"] != "experimental":
            continue
        score = agreement("13C", shifts, rec["peaks"])
        cid = rec["compound_id"]
        if cid not in best or score["score"] > best[cid][0]["score"]:
            best[cid] = (score, rec)
    ranked = sorted(best.items(), key=lambda kv: -kv[1][0]["score"])[:limit]
    return [{"smiles": isomers[cid]["smiles"],
             "inchikey": isomers[cid]["inchikey"],
             "evidence": "predicted spectrum of the query vs experimental "
                         "spectrum of a known isomer",
             "source": rec["source"], "record_id": rec["record_id"],
             "doi": rec["doi"] or None, "agreement": score}
            for cid, (score, rec) in ranked]


def _peak_list(info, nucleus):
    ref = info["reference"][nucleus]
    if ref is not None:
        return centres(nucleus, ref["peaks"]), "experimental"
    return [p["shift"] for p in info["predicted"][nucleus]], "predicted"


def _mixture(target, others, solvent, log):
    """Where signals of other species or the solvent sit on the target's."""
    rows = []
    for nucleus in NUCLEI:
        mine, mine_kind = _peak_list(target, nucleus)
        for other in others:
            if other.get("status") == "invalid_input":
                continue
            theirs, kind = _peak_list(other, nucleus)
            hits = overlaps(nucleus, mine, theirs)
            clash = {y for _, y in hits}
            rows.append({
                "nucleus": nucleus, "species": other["smiles"],
                "role": other["role"], "evidence": kind,
                "target_evidence": mine_kind,
                "overlapping": sorted({(round(x, 2), round(y, 2))
                                       for x, y in hits}, reverse=True),
                "diagnostic": sorted({round(y, 2) for y in theirs
                                      if y not in clash}, reverse=True),
            })
        if solvent:
            theirs = solvents.RESIDUAL[solvent][nucleus]
            hits = overlaps(nucleus, mine, theirs)
            rows.append({
                "nucleus": nucleus, "species": solvent, "role": "solvent",
                "evidence": "reference table", "target_evidence": mine_kind,
                "overlapping": sorted({(round(x, 2), round(y, 2))
                                       for x, y in hits}, reverse=True),
                "diagnostic": [],
            })
    return rows


def run(index, target, reactants=(), reagents=(), impurities=(), solvent=None,
        log=None):
    """Run the full check and return the report as a plain dict."""
    log = log or Log()
    reaction = None
    if ">" in target:
        reaction = target
        lhs, agents, products = parse_reaction(target)
        if not products:
            raise ValueError("the reaction has no product")
        reactants = list(lhs) + list(reactants)
        reagents = list(agents) + list(reagents)
        impurities = list(products[1:]) + list(impurities)
        target = products[0]
        if len(products) > 1:
            log("multiple_products", reaction,
                "first product taken as the target, the rest as impurities")

    solvent_key = None
    if solvent:
        solvent_key = solvents.resolve(solvent)
        if solvent_key is None:
            log("solvent_unknown", solvent,
                "no residual-peak entry for this solvent; solvent check skipped")

    for name in index.meta.get("missing_sources", []):
        log("source_not_indexed", name,
            f"{SOURCES[name].filename} was not available when the index was built")

    info, mol = describe(index, target, "target", log)
    report = {
        "tool": {"name": "nmr-novelty-check", "version": __version__},
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "query": {"target": target, "reaction": reaction,
                  "reactants": list(reactants), "reagents": list(reagents),
                  "impurities": list(impurities), "solvent": solvent},
    }

    if mol is None:
        report["verdict"] = {"status": "invalid_input", "confidence": 0.0,
                             "summary": VERDICT_TEXT["invalid_input"],
                             "reasons": []}
        report["target"] = info
        components = []
    else:
        status, conf, reasons = _verdict(index, info)
        report["verdict"] = {"status": status, "confidence": conf,
                             "summary": VERDICT_TEXT[status], "reasons": reasons}
        counts = {}
        for nucleus in NUCLEI:
            views = info["records"][nucleus]
            exp = [v for v in views
                   if v["match"] == "exact" and v["evidence"] == "experimental"]
            counts[nucleus] = {
                "exact_experimental": len(exp),
                "exact_calculated": sum(v["match"] == "exact" and
                                        v["evidence"] == "calculated"
                                        for v in views),
                "variant": sum(v["match"] == "variant" for v in views),
                "documents": len({v["doi"] for v in exp if v["doi"]}),
            }
            if len(views) > MAX_LISTED:
                log("records_truncated", f"{info['smiles']} {nucleus}",
                    f"{len(views)} records found, the {MAX_LISTED} best listed")
        report["target"] = {
            "smiles": info["smiles"], "inchikey": info["inchikey"],
            "formula": info["formula"], "atom_map": info["atom_map"],
            "record_counts": counts,
            "prediction_independent": info["prediction_independent"],
            "predicted": info["predicted"],
            "evidence": {n: info["records"][n][:MAX_LISTED] for n in NUCLEI},
            "assignment": {n: _assignment(info, mol, n) for n in NUCLEI},
        }
        for nucleus in NUCLEI:
            if info["reference"][nucleus] is None and info["predicted"][nucleus]:
                log("assignment_predicted_only", f"{info['smiles']} {nucleus}",
                    "no experimental record to assign; predicted shifts listed")
        if status in ("not_found", "listed_without_experiment"):
            report["target"]["nearest_known_isomers"] = _neighbours(index, info)

        components, full = [], []
        for role, items in (("reactant", reactants), ("reagent", reagents),
                            ("impurity", impurities)):
            for smiles in items:
                other, other_mol = describe(index, smiles, role, log)
                full.append(other)
                if other_mol is None:
                    components.append(other)
                    continue
                o_status, o_conf, _ = _verdict(index, other)
                components.append({
                    "role": role, "smiles": other["smiles"],
                    "inchikey": other["inchikey"], "formula": other["formula"],
                    "status": o_status, "confidence": o_conf,
                    "peaks": {n: {"values": _peak_list(other, n)[0],
                                  "evidence": _peak_list(other, n)[1],
                                  "source": (other["reference"][n] or {}).get("source"),
                                  "record_id": (other["reference"][n] or {}).get("record_id"),
                                  "doi": (other["reference"][n] or {}).get("doi")}
                              for n in NUCLEI},
                })
        report["components"] = components
        report["mixture"] = _mixture(info, full, solvent_key, log) \
            if (full or solvent_key) else []

    report["parameters"] = {
        "agreement_threshold": AGREE, "overlap_window_ppm": OVERLAP,
        "identity": "InChIKey (exact) / first InChIKey block (variant)",
    }
    report["data"] = {
        "index_built": index.meta["built"],
        "compounds": index.meta["compounds"],
        "sources": index.meta["sources"],
        "predictor": index.meta["predictor"],
    }
    report["degradations"] = log.entries
    report["limitations"] = LIMITATIONS
    return report
