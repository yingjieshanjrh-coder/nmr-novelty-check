"""Command line entry point."""
import argparse
import os
import re
import sys
from pathlib import Path

from . import __version__
from .check import Log, run
from .index import Index, build, index_path
from .report import to_json, to_markdown
from .sources import SOURCES, fetch_all

DEFAULT_DATA = os.environ.get("NMRCHECK_DATA", "data")


def prepare(data_dir, log, from_dir=None, offline=False, rebuild=False,
            limit=None, quiet=False):
    """Fetch what is missing and build the index if there is none."""
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    if index_path(data_dir).exists() and not rebuild:
        return Index(data_dir)
    files = fetch_all(data_dir, from_dir=from_dir, offline=offline, log=log)
    if not files:
        raise SystemExit(
            "no data source could be obtained; check the network or pass "
            "--from-dir with a folder that holds the downloaded files")
    if not quiet:
        print("building index (first run only; the full NMRexp set takes a "
              "while)", file=sys.stderr)
    build(data_dir, files, log=log, limit=limit, quiet=quiet)
    return Index(data_dir)


def slug(text):
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")[:60] or "report"


def write_report(report, out_dir, name):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for ext, render in (("json", to_json), ("md", to_markdown)):
        path = out_dir / f"{name}.{ext}"
        path.write_text(render(report), encoding="utf-8")
        paths[ext] = path
    return paths


def add_data_options(parser):
    parser.add_argument("--data-dir", default=DEFAULT_DATA,
                        help="where raw files and the index live (default: ./data)")
    parser.add_argument("--from-dir",
                        help="folder with already downloaded source files")
    parser.add_argument("--offline", action="store_true",
                        help="never touch the network")


def cmd_fetch(args):
    log = Log()
    files = fetch_all(args.data_dir, from_dir=args.from_dir, offline=args.offline,
                      log=log)
    for name in SOURCES:
        state = files.get(name, "missing")
        print(f"{name:12s} {state}")
    for entry in log.entries:
        print(f"  note: {entry['subject']}: {entry['detail']}", file=sys.stderr)
    return 0 if files else 1


def cmd_build(args):
    log = Log()
    index = prepare(args.data_dir, log, from_dir=args.from_dir,
                    offline=args.offline, rebuild=True, limit=args.limit)
    for name, info in index.meta["sources"].items():
        print(f"{name:12s} {info['spectra']} spectra")
    print(f"compounds    {index.meta['compounds']}")
    for entry in log.entries:
        print(f"  note: {entry['subject']}: {entry['detail']}", file=sys.stderr)
    return 0


def cmd_check(args):
    log = Log()
    index = prepare(args.data_dir, log, from_dir=args.from_dir,
                    offline=args.offline)
    report = run(index, args.target, reactants=args.reactant,
                 reagents=args.reagent, impurities=args.impurity,
                 solvent=args.solvent, log=log)
    paths = write_report(report, args.out, args.name or slug(args.target))
    v = report["verdict"]
    print(f"{v['status']} (confidence {v['confidence']:.2f})")
    print(v["summary"])
    for kind, path in paths.items():
        print(f"{kind}: {path}")
    return 2 if v["status"] == "invalid_input" else 0


def cmd_serve(args):
    from .web import create_app
    log = Log()
    index = prepare(args.data_dir, log, from_dir=args.from_dir,
                    offline=args.offline)
    app = create_app(index, startup_log=log)
    app.run(host=args.host, port=args.port)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="nmrcheck",
        description="Check whether NMR data for a compound are already public.")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("check", help="run a check (fetches and indexes on first use)")
    p.add_argument("target", help="SMILES of the target, or a reaction SMILES "
                                  "reactants>agents>products")
    p.add_argument("-r", "--reactant", action="append", default=[], metavar="SMILES")
    p.add_argument("-g", "--reagent", action="append", default=[], metavar="SMILES")
    p.add_argument("-i", "--impurity", action="append", default=[], metavar="SMILES")
    p.add_argument("-s", "--solvent", help="NMR solvent, e.g. CDCl3")
    p.add_argument("-o", "--out", default="reports", help="output folder")
    p.add_argument("--name", help="file stem for the report")
    add_data_options(p)
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("fetch", help="download the raw source files")
    add_data_options(p)
    p.set_defaults(func=cmd_fetch)

    p = sub.add_parser("build", help="(re)build the local index")
    p.add_argument("--limit", type=int,
                   help="index only the first N NMRexp records (for a quick trial)")
    add_data_options(p)
    p.set_defaults(func=cmd_build)

    p = sub.add_parser("serve", help="start the web form")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8000)
    add_data_options(p)
    p.set_defaults(func=cmd_serve)

    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
