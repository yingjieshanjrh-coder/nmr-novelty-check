"""Render a report dict as JSON, Markdown or HTML."""
import html
import json
import re


def to_json(report):
    return json.dumps(report, indent=2, ensure_ascii=False)


def _fmt(value, digits=2):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def _peaks(nucleus, peaks, limit=40):
    if nucleus == "1H":
        parts = []
        for centre, low, high, mult, n_h in peaks[:limit]:
            pos = f"{high:.2f}-{low:.2f}" if high != low else f"{centre:.2f}"
            extra = ", ".join(x for x in (mult, f"{n_h}H" if n_h else "") if x)
            parts.append(f"{pos} ({extra})" if extra else pos)
    else:
        parts = [f"{p:.2f}" for p in peaks[:limit]]
    more = f" ... (+{len(peaks) - limit})" if len(peaks) > limit else ""
    return ", ".join(parts) + more


def _table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(" --- " for _ in header) + "|"]
    for row in rows:
        out.append("| " + " | ".join(str(c).replace("|", "\\|") for c in row) + " |")
    return "\n".join(out)


def _ref(view):
    if view.get("doi"):
        return f"[{view['doi']}](https://doi.org/{view['doi']})"
    if view.get("source") == "nmrshiftdb2":
        return f"nmrshiftdb2 molecule {view['record_id'].split('/')[0]}"
    return f"{view.get('source')} #{view.get('record_id')}"


def to_markdown(report):
    v = report["verdict"]
    q = report["query"]
    lines = ["# NMR novelty check", ""]
    lines.append(f"Query: `{q['reaction'] or q['target']}`  ")
    lines.append(f"Generated: {report['generated']} "
                 f"(nmr-novelty-check {report['tool']['version']})")
    lines += ["", "## Verdict", "",
              f"**{v['status']}**, confidence {v['confidence']:.2f}", "",
              v["summary"], ""]
    for reason in v["reasons"]:
        lines.append(f"- {reason}")

    t = report.get("target", {})
    if "inchikey" in t:
        lines += ["", "## Target", "",
                  _table(["field", "value"], [
                      ["canonical SMILES", f"`{t['smiles']}`"],
                      ["InChIKey", t["inchikey"]],
                      ["formula", t["formula"]],
                      ["atom numbering", f"`{t['atom_map']}`"],
                  ])]
        rows = [[n, c["exact_experimental"], c["documents"],
                 c["exact_calculated"], c["variant"]]
                for n, c in t["record_counts"].items()]
        lines += ["", "Records found:", "",
                  _table(["nucleus", "exact, experimental", "distinct DOIs",
                          "exact, calculated", "stereo/isotope variants"], rows)]

        for nucleus, views in t["evidence"].items():
            if not views:
                continue
            lines += ["", f"### Evidence, {nucleus}", ""]
            rows = [[x["evidence"], x["match"], x["source"], _ref(x),
                     x["solvent"] or "-", _fmt(x["agreement"]["score"]),
                     _fmt(x["agreement"]["mae"]), _peaks(nucleus, x["peaks"])]
                    for x in views]
            lines.append(_table(["evidence", "match", "source", "reference",
                                 "solvent", "agreement", "MAE / ppm", "peaks / ppm"],
                                rows))

        for nucleus, rows in t["assignment"].items():
            if not rows:
                continue
            lines += ["", f"### Assignment, {nucleus}", ""]
            body = [[r["atoms"] or "unassigned", _fmt(r["observed"]),
                     _fmt(r["predicted"]), _fmt(r["delta"]), _fmt(r["count"]),
                     _fmt(r["radius"]), _fmt(r["support"])] for r in rows]
            lines.append(_table(["atoms", "observed", "predicted", "delta",
                                 "nuclei", "env. radius", "training hits"], body))
        if not t["prediction_independent"]:
            lines += ["", "The predictor was trained on data that include this "
                          "compound; observed and predicted values are not "
                          "independent."]

        near = t.get("nearest_known_isomers")
        if near is not None:
            lines += ["", "### Closest known isomers (13C)", ""]
            if near:
                lines.append(_table(
                    ["SMILES", "agreement", "MAE / ppm", "reference"],
                    [[f"`{n['smiles']}`", _fmt(n["agreement"]["score"]),
                      _fmt(n["agreement"]["mae"]), _ref(n)] for n in near]))
            else:
                lines.append("No isomer with an experimental 13C spectrum is indexed.")

    comps = report.get("components") or []
    if comps:
        lines += ["", "## Other species", ""]
        rows = []
        for c in comps:
            if c.get("status") == "invalid_input":
                rows.append([c["role"], f"`{c['input']}`", "invalid_input", "-",
                             "-", "-"])
                continue
            rows.append([c["role"], f"`{c['smiles']}`", c["status"],
                         _fmt(c["confidence"]),
                         f"{c['peaks']['13C']['evidence']}: "
                         + (_peaks("13C", c["peaks"]["13C"]["values"]) or "-"),
                         f"{c['peaks']['1H']['evidence']}: "
                         + (_peaks("13C", c["peaks"]["1H"]["values"]) or "-")])
        lines.append(_table(["role", "SMILES", "status", "confidence",
                             "13C / ppm", "1H / ppm"], rows))

    mix = report.get("mixture") or []
    if mix:
        lines += ["", "## Signal overlap with the target", ""]
        rows = []
        for m in mix:
            over = ", ".join(f"{x:.2f}~{y:.2f}" for x, y in m["overlapping"]) or "none"
            diag = ", ".join(f"{x:.2f}" for x in m["diagnostic"][:12]) or "-"
            rows.append([m["nucleus"], f"`{m['species']}`", m["role"],
                         m["evidence"], over, diag])
        lines.append(_table(["nucleus", "species", "role", "evidence",
                             "overlap (target~species)", "signals clear of target"],
                            rows))

    lines += ["", "## Data sources", ""]
    d = report["data"]
    rows = [[name, s["file"], s["spectra"], s["license"], s["md5"][:12],
             s["citation"]] for name, s in d["sources"].items()]
    lines.append(_table(["source", "file", "spectra indexed", "licence", "md5",
                         "citation"], rows))
    lines += ["", f"Index built {d['index_built']}, {d['compounds']} compounds."]
    hold = d["predictor"].get("holdout", {})
    if hold:
        bits = [f"{n}: MAE {s['all']['mae']} ppm over {s['all']['n']} signals"
                for n, s in hold.items()]
        lines.append("Predictor hold-out error (10% of nmrshiftdb2 compounds): "
                     + "; ".join(bits) + ".")

    lines += ["", "## Degradation log", ""]
    if report["degradations"]:
        lines.append(_table(["code", "subject", "detail"],
                            [[e["code"], f"`{e['subject']}`", e["detail"]]
                             for e in report["degradations"]]))
    else:
        lines.append("Nothing was degraded in this run.")

    lines += ["", "## Limitations", ""]
    lines += [f"- {x}" for x in report["limitations"]]
    return "\n".join(lines) + "\n"


def _inline(text):
    text = html.escape(text)
    text = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r'<a href="\2">\1</a>', text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", text)


def to_html(report):
    """Small Markdown-subset renderer; enough for the report layout above."""
    out, rows, items = [], [], []

    def flush():
        nonlocal rows, items
        if rows:
            head, body = rows[0], rows[2:]
            out.append("<table><thead><tr>" + "".join(
                f"<th>{_inline(c)}</th>" for c in head) + "</tr></thead><tbody>")
            for r in body:
                out.append("<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in r)
                           + "</tr>")
            out.append("</tbody></table>")
            rows = []
        if items:
            out.append("<ul>" + "".join(f"<li>{_inline(i)}</li>" for i in items)
                       + "</ul>")
            items = []

    for line in to_markdown(report).splitlines():
        if line.startswith("|"):
            cells = [c.strip().replace("\\|", "|")
                     for c in re.split(r"(?<!\\)\|", line.strip()[1:-1])]
            rows.append(cells)
            continue
        if line.startswith("- "):
            items.append(line[2:])
            continue
        flush()
        m = re.match(r"(#+) (.*)", line)
        if m:
            level = len(m.group(1))
            out.append(f"<h{level}>{_inline(m.group(2))}</h{level}>")
        elif line.strip():
            out.append(f"<p>{_inline(line.strip())}</p>")
    flush()
    return "\n".join(out)
