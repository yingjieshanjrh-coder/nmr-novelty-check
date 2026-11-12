"""From the two input files to a classified peak list."""
from . import impurities
from .assign import Group, Settings, Species, analyse
from .predictor import build_species


def trace_species(name, solvent):
    """A declared lab solvent as a species, or as loose lines if it has no
    fixed proton ratio (grease, TMS)."""
    lines = impurities.trace_lines(name, solvent)
    if any(protons is None for _, _, protons in lines):
        return None, [{"name": name, "kind": "trace", "ppm": ppm,
                       "window": impurities.TRACE_WINDOW}
                      for ppm, _, _ in lines]
    groups = [Group(label=f"{mult}, {protons}H", shift=ppm, protons=protons,
                    window=impurities.TRACE_WINDOW,
                    tight=impurities.TRACE_WINDOW)
              for ppm, mult, protons in lines]
    return Species(key=f"trace:{name}", role="impurity", name=name,
                   groups=groups), []


def run(targets, peaks, index=None, cfg=None, tol=0.08, predict=build_species):
    """Classify `peaks` for the species in `targets`.

    `predict` can be replaced in tests to avoid needing RDKit and the index.
    """
    cfg = cfg or Settings()
    product = predict(targets.product.smiles, "product", "product",
                      targets.product.name, index, tol)
    reactants = [predict(r.smiles, f"reactant{i}", "reactant", r.name,
                         index, tol)
                 for i, r in enumerate(targets.reactants, start=1)]

    lines = impurities.solvent_lines(targets.solvent)
    extra = []
    for name in targets.extra_solvents:
        species, loose = trace_species(name, targets.solvent)
        if species:
            extra.append(species)
        lines += loose
    for item in targets.extra_peaks:
        lines.append({"name": item["name"], "kind": "trace",
                      "ppm": item["ppm"], "window": impurities.TRACE_WINDOW})

    result = analyse(peaks, product, reactants, lines=lines, impurities=extra,
                     cfg=cfg, reference=targets.reference,
                     solvent=targets.solvent)
    if index is None:
        result.warnings.append(
            "没有找到父项目的索引（data/index.sqlite），预测位移退化为按原子类型给的"
            "固定值，窗口很宽，归属仅供参考。请先在父文件夹运行 nmrcheck build。")
    weak = [g for s in [product] + reactants for g in s.groups
            if g.radius < 2 and not g.labile]
    if weak and index is not None:
        result.warnings.append(
            f"有 {len(weak)} 组信号的预测只基于很小的原子环境（半径 < 2 根键），"
            "搜索窗口相应放宽了，这些归属要多留个心眼。")
    return result
