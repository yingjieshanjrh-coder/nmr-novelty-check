"""Command line entry point."""
import argparse
import sys
from pathlib import Path

from . import __version__
from .assign import Settings
from .inputs import InputError, load_peaks, load_targets
from .parent import open_index
from .pipeline import run
from .report import to_json, to_markdown, to_text


def parser():
    p = argparse.ArgumentParser(
        prog="peakassign",
        description="Sort the peaks of a crude 1H NMR spectrum into product, "
                    "leftover reactants, solvent and unexplained signals.")
    p.add_argument("targets", help="targets.yaml: product, reactants, solvent")
    p.add_argument("peaks", help="measured_peaks.csv: ppm, integral[, multiplicity]")
    p.add_argument("-o", "--output", help="Markdown report (default: "
                   "report.md next to the peak list)")
    p.add_argument("--json", help="also write the result as JSON")
    p.add_argument("--data-dir", help="folder holding the parent project's "
                   "index.sqlite (default: ../data)")
    p.add_argument("--tol", type=float, default=0.08,
                   help="ppm tolerance for a confident match (default 0.08); "
                        "the search window adds the predictor's uncertainty")
    p.add_argument("--unknown-threshold", type=float, default=0.5,
                   help="unexplained integral, in protons, that raises the "
                        "unknown warning (default 0.5)")
    p.add_argument("--residual-alert", type=float, default=0.20,
                   help="reactant:product molar ratio flagged as a lot "
                        "(default 0.20)")
    p.add_argument("--strict-overlap", action="store_true",
                   help="list every product peak that also lies in a reactant "
                        "window as ambiguous, even if the integral fits")
    p.add_argument("--version", action="version", version=__version__)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    cfg = Settings(unknown_h=args.unknown_threshold,
                   residual_alert=args.residual_alert,
                   strict_overlap=args.strict_overlap)
    try:
        targets = load_targets(args.targets)
        peaks = load_peaks(args.peaks)
        index = open_index(args.data_dir)
        try:
            result = run(targets, peaks, index=index, cfg=cfg, tol=args.tol)
        finally:
            if index is not None:
                index.close()
    except ImportError as exc:
        print(f"error: {exc}\nRDKit and the parent nmrcheck package are "
              "needed; run `pip install -e .` in the parent folder first.",
              file=sys.stderr)
        return 2
    except (InputError, ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    output = Path(args.output) if args.output else \
        Path(args.peaks).with_name("report.md")
    output.write_text(to_markdown(result), encoding="utf-8")
    if args.json:
        Path(args.json).write_text(to_json(result), encoding="utf-8")
    print(to_text(result))
    print(f"report: {output}")
    return 0
