"""Local SQLite index over all sources."""
import json
import multiprocessing as mp
import os
import sqlite3
import sys
import time
import zlib
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .chem import identify, identify_smiles
from .parsers import READERS
from .predict import MAX_RADIUS, ShiftTable, predict
from .sources import SOURCES, file_md5

SCHEMA = """
CREATE TABLE compounds (
    id INTEGER PRIMARY KEY,
    inchikey TEXT NOT NULL UNIQUE,
    skeleton TEXT NOT NULL,
    smiles TEXT NOT NULL,
    formula TEXT NOT NULL
);
CREATE TABLE spectra (
    id INTEGER PRIMARY KEY,
    compound_id INTEGER NOT NULL,
    source TEXT NOT NULL,
    record_id TEXT NOT NULL,
    nucleus TEXT NOT NULL,
    origin TEXT NOT NULL,
    solvent TEXT,
    frequency TEXT,
    doi TEXT,
    peaks TEXT NOT NULL
);
CREATE TABLE shift_env (
    nucleus TEXT NOT NULL,
    radius INTEGER NOT NULL,
    env INTEGER NOT NULL,
    n INTEGER NOT NULL,
    shift REAL NOT NULL,
    spread REAL NOT NULL,
    PRIMARY KEY (nucleus, radius, env)
) WITHOUT ROWID;
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
"""

INDEXES = """
CREATE INDEX compounds_skeleton ON compounds (skeleton);
CREATE INDEX compounds_formula ON compounds (formula);
CREATE INDEX spectra_compound ON spectra (compound_id);
"""

HOLDOUT_EVERY = 10


def index_path(data_dir):
    return Path(data_dir) / "index.sqlite"


class Builder:
    def __init__(self, db):
        self.db = db
        self.compound_ids = {}
        self.counts = {}

    def compound(self, inchikey, smiles, formula):
        cid = self.compound_ids.get(inchikey)
        if cid is None:
            cid = len(self.compound_ids) + 1
            self.compound_ids[inchikey] = cid
            self.db.execute("INSERT INTO compounds VALUES (?,?,?,?,?)",
                            (cid, inchikey, inchikey[:14], smiles, formula))
        return cid

    def spectrum(self, cid, source, rec):
        self.db.execute(
            "INSERT INTO spectra (compound_id, source, record_id, nucleus, origin,"
            " solvent, frequency, doi, peaks) VALUES (?,?,?,?,?,?,?,?,?)",
            (cid, source, rec["record_id"], rec["nucleus"], rec["origin"],
             rec["solvent"], rec["frequency"], rec["doi"],
             json.dumps(rec["peaks"], separators=(",", ":"))))
        self.counts[source] = self.counts.get(source, 0) + 1


def _is_holdout(inchikey):
    return zlib.crc32(inchikey.encode()) % HOLDOUT_EVERY == 0


def _evaluate(table, held):
    """Mean absolute error of the predictor on compounds it has not seen."""
    errors = {}
    for mol, nucleus, assigned in held:
        truth = {}
        for atom, shift in assigned:
            truth.setdefault(atom, []).append(shift)
        for peak in predict(mol, nucleus, table.lookup):
            obs = [s for a in peak["atoms"] for s in truth.get(a, [])]
            if not obs:
                continue
            err = abs(peak["shift"] - sum(obs) / len(obs))
            errors.setdefault((nucleus, peak["radius"]), []).append(err)
    stats = {}
    for (nucleus, radius), errs in errors.items():
        stats.setdefault(nucleus, {})[str(radius)] = {
            "n": len(errs), "mae": round(sum(errs) / len(errs), 3)}
    for nucleus in list(stats):
        stats[nucleus] = dict(sorted(stats[nucleus].items(), reverse=True))
        errs = [e for (nuc, _), v in errors.items() if nuc == nucleus for e in v]
        stats[nucleus]["all"] = {"n": len(errs),
                                 "mae": round(sum(errs) / len(errs), 3)}
    return stats


def _load_nmrshiftdb2(builder, path, progress):
    table = ShiftTable()
    held = []
    seen = 0
    for rec in READERS["nmrshiftdb2"](path):
        mol = rec.pop("mol")
        assigned = rec.pop("assigned")
        ident = identify(mol)
        if ident is None:
            continue
        cid = builder.compound(ident.inchikey, ident.smiles, ident.formula)
        builder.spectrum(cid, "nmrshiftdb2", rec)
        if rec["origin"] == "experimental" and rec["nucleus"] in ("13C", "1H"):
            if _is_holdout(ident.inchikey):
                held.append((mol, rec["nucleus"], assigned))
            else:
                table.add(mol, rec["nucleus"], assigned)
        seen += 1
        if seen % 10000 == 0:
            progress(f"nmrshiftdb2: {seen} spectra")
    stats = _evaluate(table, held)
    for mol, nucleus, assigned in held:
        table.add(mol, nucleus, assigned)
    builder.db.executemany("INSERT INTO shift_env VALUES (?,?,?,?,?,?)",
                           table.rows())
    return stats


def _load_smiles_source(builder, name, path, progress, workers, limit=None):
    cache = {}
    pool = mp.get_context("spawn").Pool(workers) if workers > 1 else None
    seen = 0
    try:
        pending = []

        def flush():
            fresh = list({r["smiles"] for r in pending if r["smiles"] not in cache})
            if fresh:
                if pool:
                    idents = pool.map(identify_smiles, fresh, chunksize=2000)
                else:
                    idents = [identify_smiles(s) for s in fresh]
                for smiles, ident in zip(fresh, idents):
                    cache[smiles] = (builder.compound(*ident) if ident else None)
            for rec in pending:
                cid = cache[rec["smiles"]]
                if cid is not None:
                    builder.spectrum(cid, name, rec)
            pending.clear()

        for rec in READERS[name](path):
            pending.append(rec)
            seen += 1
            if len(pending) >= 100000:
                flush()
                progress(f"{name}: {seen} spectra")
            if limit and seen >= limit:
                break
        flush()
    finally:
        if pool:
            pool.close()
            pool.join()
    return sum(1 for v in cache.values() if v is None)


def build(data_dir, files, log=None, workers=None, limit=None, quiet=False):
    """Build the index from {source name: raw file path}."""
    note = log or (lambda *a, **k: None)
    target = index_path(data_dir)
    tmp = target.with_suffix(".building")
    if tmp.exists():
        tmp.unlink()
    workers = workers or max(1, (os.cpu_count() or 2))
    started = time.time()

    def progress(msg):
        if not quiet:
            print(f"  [{time.time() - started:6.0f}s] {msg}", file=sys.stderr)

    db = sqlite3.connect(tmp)
    db.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;")
    db.executescript(SCHEMA)
    builder = Builder(db)
    meta = {"version": __version__,
            "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "sources": {}}

    predictor_stats = {}
    if "nmrshiftdb2" in files:
        predictor_stats = _load_nmrshiftdb2(builder, files["nmrshiftdb2"], progress)
    else:
        note("predictor_untrained", "nmrshiftdb2",
             "no assigned spectra available; predictions use fixed atom-type values")
    for name, path in files.items():
        if name == "nmrshiftdb2":
            continue
        try:
            bad = _load_smiles_source(builder, name, path, progress, workers, limit)
        except Exception as exc:
            note("source_unreadable", name, f"{type(exc).__name__}: {exc}")
            continue
        if bad:
            note("structures_skipped", name,
                 f"{bad} distinct SMILES could not be parsed and were left out")
        if limit:
            note("source_truncated", name, f"only the first {limit} records indexed")

    progress("creating indexes")
    db.executescript(INDEXES)
    for name, path in files.items():
        src = SOURCES[name]
        stamp = Path(str(path) + ".md5")
        md5 = stamp.read_text().strip() if stamp.exists() else file_md5(path)
        meta["sources"][name] = {
            "file": src.filename, "md5": md5, "bytes": os.path.getsize(path),
            "spectra": builder.counts.get(name, 0), "license": src.license,
            "citation": src.citation, "homepage": src.homepage,
        }
    meta["missing_sources"] = [n for n in SOURCES if n not in files]
    meta["compounds"] = len(builder.compound_ids)
    meta["predictor"] = {"max_radius": MAX_RADIUS, "holdout": predictor_stats}
    db.execute("INSERT INTO meta VALUES ('info', ?)", (json.dumps(meta),))
    db.commit()
    db.close()
    tmp.replace(target)
    progress(f"index ready: {meta['compounds']} compounds")
    return meta


class Index:
    def __init__(self, data_dir):
        self.path = index_path(data_dir)
        self.db = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True,
                                  check_same_thread=False)
        self.db.row_factory = sqlite3.Row
        row = self.db.execute("SELECT value FROM meta WHERE key='info'").fetchone()
        self.meta = json.loads(row["value"])

    def close(self):
        self.db.close()

    def lookup_env(self, nucleus, radius, ident):
        row = self.db.execute(
            "SELECT n, shift FROM shift_env WHERE nucleus=? AND radius=? AND env=?",
            (nucleus, radius, ident)).fetchone()
        return (row["n"], row["shift"]) if row else None

    def compounds(self, inchikey=None, skeleton=None, formula=None, limit=2000):
        if inchikey:
            sql, arg = "inchikey=?", inchikey
        elif skeleton:
            sql, arg = "skeleton=?", skeleton
        else:
            sql, arg = "formula=?", formula
        return self.db.execute(
            f"SELECT * FROM compounds WHERE {sql} LIMIT ?", (arg, limit)).fetchall()

    def spectra(self, compound_ids, nucleus=None):
        if not compound_ids:
            return []
        marks = ",".join("?" * len(compound_ids))
        sql = f"SELECT * FROM spectra WHERE compound_id IN ({marks})"
        args = list(compound_ids)
        if nucleus:
            sql += " AND nucleus=?"
            args.append(nucleus)
        out = []
        for row in self.db.execute(sql, args):
            item = dict(row)
            item["peaks"] = json.loads(item["peaks"])
            out.append(item)
        return out
