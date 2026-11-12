"""Report rendering: Markdown for reading, JSON for scripts."""
import json

from . import __version__, impurities

ALERT = "发现疑似未知副产物信号，建议重点排查"
STATUS_NOTE = {
    "excess": "积分比预测多 {left:.2f}H，下面可能还压着别的信号",
    "split": "与相邻的峰合并后积分才对得上（同一组氢被分成了几个峰）",
    "labile": "活泼氢，位移和积分随浓度、含水量变化，仅供参考",
}
SIMPLE = {"s", "d", "t", "q", "sep"}
SCALE_NOTE = {
    "product": "由产物信号整体拟合得到",
    "reference": "按 targets.yaml 中指定的参考峰",
    "total": "产物未匹配上，按总面积粗略估计",
}


def _signal(part):
    labels = " + ".join(part.labels)
    shifts = " / ".join(f"{s:.2f}" for s in part.predicted)
    return labels, shifts


def _notes(state, part):
    notes = []
    if part.status in STATUS_NOTE:
        notes.append(STATUS_NOTE[part.status].format(left=state.left_h))
    if state.also:
        notes.append("也在 " + "、".join(state.also) + " 的预测窗口内，"
                     "但积分没有多出来的面积")
    return "；".join(notes)


def _table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(" --- " for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return out


def _assigned_rows(states, key):
    rows = []
    for state in states:
        for part in state.parts:
            if part.species != key:
                continue
            labels, shifts = _signal(part)
            rows.append([f"{state.ppm:.2f}", state.multiplicity or "-",
                         f"{part.h:.2f}", labels, shifts, f"{part.protons:g}",
                         f"{part.delta:.2f}", part.confidence,
                         _notes(state, part) or "-"])
    return rows


ASSIGNED_HEADER = ["实测 ppm", "峰形", "积分 (H)", "归属建议", "预测 ppm",
                   "预测 H", "偏差", "置信度", "备注"]


def _unknown_hints(state, solvent):
    hints = list(state.hints)
    seen = state.multiplicity.strip().lower()
    near = []
    for name, shift, mult, lines in impurities.nearby_trace(solvent, state.ppm):
        if seen in SIMPLE and mult in SIMPLE and seen != mult:
            continue
        more = f"，共 {lines} 组峰" if lines > 1 else ""
        near.append(f"{name}（{shift:.2f}, {mult}{more}）")
    if near:
        hints.append(f"位移接近 {solvent} 中的 " + "、".join(near)
                     + "；如果实验里用过，在 targets.yaml 的 extra_solvents "
                     "里声明后重跑")
    return "；".join(hints) or "不在任何候选物的预测范围内"


def _ambiguous_text(state):
    if state.conflict:
        return "冲突", state.conflict
    bits = []
    for part in state.parts:
        labels, shifts = _signal(part)
        if part.role == "solvent":
            bits.append(f"{part.name}（{shifts}）约 {part.h:.2f}H")
        else:
            bits.append(f"{part.name} 的 {labels}（预测 {shifts}）约 {part.h:.2f}H")
    for item in state.overlaps:
        bits.append(f"已检出的 {item['name']} 的 {item['label']}（预测 "
                    f"{item['shift']:.2f}）按比例应在这里贡献约 {item['h']:.2f}H")
    if not state.overlaps and len({p.species for p in state.parts}) < 2:
        bits += [f"也在 {text} 的窗口内" for text in state.also]
    text = "；".join(bits)
    if state.left_h:
        text += f"；另有 {state.left_h:.2f}H 无法解释"
    return "重叠", text


def _species_title(report):
    species = report.species
    code = f"（`{species.smiles}`）" if species.smiles and \
        species.smiles != species.name else ""
    return f"{species.name}{code}"


def _group_list(items):
    return "、".join(f"{g.label} {g.shift:.2f} ppm" for g in items)


def _reactant_block(report, result):
    lines = [f"### {_species_title(report)}", ""]
    groups = _group_list
    if report.status == "detected":
        level = f"检出，约为产物的 **{report.ratio * 100:.0f} mol%**，" \
                f"占全谱积分的 {report.share * 100:.0f}%。"
        if report.alert:
            level += " **原料残留较多。**"
        lines += [level, ""]
        lines += _table(ASSIGNED_HEADER,
                        _assigned_rows(result.peaks, report.species.key))
        if report.missing:
            lines += ["", f"没有找到的信号：{groups(report.missing)}。"
                          "归属的把握因此打了折扣。"]
        if report.masked:
            lines += ["", f"被其他峰挡住的信号：{groups(report.masked)}，"
                          "见第 5 节。"]
    elif report.status == "doubtful":
        lines += ["有个别峰落在它的预测窗口内，但多数信号"
                  f"（{groups(report.missing)}）没有出现，**未按该原料归属**，"
                  "相关的峰留在第 4 节。"]
    elif report.status == "masked":
        lines += [f"未检出。它的信号（{groups(report.masked)}）预测位置都被已"
                  "归属的峰占着，单凭位移排除不了；不过这些峰的积分没有多出来的"
                  "面积，即使有残留量也很少。"]
    else:
        text = "未检出。"
        if report.missing:
            text += f"预测信号 {groups(report.missing)} 附近没有对得上的峰。"
        lines += [text]
    return lines + [""]


def to_markdown(result, title="NMR 粗产物谱峰归属报告"):
    cfg = result.settings
    product = result.product
    unknown = result.category("unknown")
    counts = {name: len(result.category(name)) for name in
              ("product", "reactant", "solvent", "unknown", "ambiguous",
               "minor")}

    out = [f"# {title}", ""]
    if result.unknown_alert:
        where = [f"{s.ppm:.2f} ppm（{s.left_h:.2f}H）"
                 for s in unknown + result.hidden]
        out += [f"> ⚠️ **{ALERT}：{'、'.join(where)}**", ""]
    else:
        out += [f"> 没有发现积分超过 {cfg.unknown_h:g}H 的未知信号。", ""]
    for text in result.warnings:
        out += [f"> 注意：{text}", ""]

    reactants = "、".join(_species_title(r) for r in result.reactants) or "无"
    out += _table(["项目", "内容"], [
        ["预期产物", _species_title(product)],
        ["反应物", reactants],
        ["氘代溶剂", result.solvent or "-"],
        ["积分归一化", f"1H = 原始积分 {result.unit:.3f}"
                      f"（{SCALE_NOTE[result.scale_source]}）"],
        ["实测峰", f"{len(result.peaks)} 个：产物 {counts['product']}，"
                  f"原料 {counts['reactant']}，溶剂/杂质 {counts['solvent']}，"
                  f"未知 {counts['unknown']}，冲突/重叠 {counts['ambiguous']}，"
                  f"未归属小峰 {counts['minor']}"],
    ])

    out += ["", "## 1. 产物特征峰（Product Assigned）", ""]
    rows = _assigned_rows(result.category("product"), product.species.key)
    if rows:
        out += _table(ASSIGNED_HEADER, rows)
    else:
        out += ["没有峰被明确归属给产物。"]
    if product.missing:
        lost = "、".join(f"{g.label}（预测 {g.shift:.2f} ppm，{g.protons}H）"
                        for g in product.missing)
        out += ["", f"**没有找到的产物信号：{lost}。** 要么产物结构不对，要么"
                    "它被别的峰盖住了（看第 5 节）。"]
    if product.species.atom_map:
        out += ["", f"原子编号：`{product.species.atom_map}`"]

    out += ["", "## 2. 未反应完的原料（Reactant Residual）", ""]
    if not result.reactants:
        out += ["没有提供反应物。", ""]
    for report in result.reactants:
        out += _reactant_block(report, result)

    out += ["## 3. 溶剂与已知杂质（Solvent / Water）", ""]
    rows = []
    for state in result.category("solvent"):
        for part in state.parts:
            _, shifts = _signal(part)
            rows.append([f"{state.ppm:.2f}", state.multiplicity or "-",
                         f"{part.h:.2f}", part.name, shifts])
    if rows:
        out += _table(["实测 ppm", "峰形", "积分 (H)", "归属", "参考 ppm"], rows)
    else:
        out += ["没有识别到溶剂峰或水峰。"]
    for report in result.impurities:
        if report.status != "detected":
            out += ["", f"声明的 {report.species.name} 没有完整出现"
                        "（各组峰需同时存在且积分成比例）。"]

    out += ["", "## 4. ⚠️ 未知副产物预警（Unassigned / Unknowns）", ""]
    if result.unknown_alert:
        out += [f"> **{ALERT}。** 下面的信号积分显著（≥ {cfg.unknown_h:g}H），"
                "产物、原料、溶剂都解释不了。", ""]
        if unknown:
            out += _table(
                ["实测 ppm", "峰形", "积分 (H)", "原始积分", "线索"],
                [[f"**{s.ppm:.2f}**", s.multiplicity or "-", f"**{s.h:.2f}**",
                  f"{s.integral:g}", _unknown_hints(s, result.solvent)]
                 for s in unknown])
        if result.hidden:
            out += ["", "已归属的峰下面多出来的面积：", ""]
            out += _table(
                ["实测 ppm", "总积分 (H)", "已解释 (H)", "多出 (H)"],
                [[f"{s.ppm:.2f}", f"{s.h:.2f}", f"{s.h - s.left_h:.2f}",
                  f"**{s.left_h:.2f}**"] for s in result.hidden])
    else:
        out += ["没有。所有积分显著的峰都有了去处。"]
    minor = result.category("minor")
    if minor:
        out += ["", f"积分不足 {cfg.unknown_h:g}H 的未归属小峰（不触发预警）："
                + "、".join(f"{s.ppm:.2f} ppm（{s.h:.2f}H）" for s in minor)
                + "。"]

    out += ["", "## 5. 冲突与重叠（Ambiguous）", ""]
    ambiguous = result.category("ambiguous")
    if ambiguous:
        rows = []
        for state in ambiguous:
            kind, text = _ambiguous_text(state)
            rows.append([f"{state.ppm:.2f}", state.multiplicity or "-",
                         f"{state.h:.2f}", kind, text])
        out += _table(["实测 ppm", "峰形", "积分 (H)", "类型", "说明"], rows)
    else:
        out += ["没有。"]

    out += ["", "## 附：预测信号", ""]
    rows = []
    for report in [product] + result.reactants:
        hit = {id(g) for g in report.matched}
        for group in report.species.groups:
            state = "已归属" if id(group) in hit else \
                ("活泼氢，未找到" if group.labile else "未归属")
            window = max(group.window, cfg.labile_window) if group.labile \
                else group.window
            rows.append([report.species.name, group.label,
                         f"{group.shift:.2f}", group.protons,
                         f"±{window:.2f}", state])
    out += _table(["物质", "信号", "预测 ppm", "H", "搜索窗口", "结果"], rows)
    out += ["", f"<sub>peakassign {__version__} · 预测位移是估算值，归属是建议，"
                "不能代替二维谱。</sub>", ""]
    return "\n".join(out)


def to_dict(result):
    def group(g):
        return {"label": g.label, "shift": g.shift, "protons": g.protons,
                "window": g.window, "labile": g.labile}

    def species(report):
        return {
            "name": report.species.name, "smiles": report.species.smiles,
            "status": report.status, "molar_ratio": report.ratio,
            "integral_h": round(report.h, 3),
            "integral_share": round(report.share, 4), "alert": report.alert,
            "matched": [group(g) for g in report.matched],
            "missing": [group(g) for g in report.missing],
            "masked": [group(g) for g in report.masked],
        }

    peaks = []
    for s in result.peaks:
        peaks.append({
            "ppm": s.ppm, "integral": s.integral, "h": round(s.h, 3),
            "multiplicity": s.multiplicity, "category": s.category,
            "unexplained_h": round(s.left_h, 3),
            "assignments": [{
                "species": p.name, "role": p.role, "signals": p.labels,
                "predicted_ppm": p.predicted, "predicted_h": p.protons,
                "integral_h": round(p.h, 3), "status": p.status,
                "delta": round(p.delta, 3), "confidence": p.confidence,
            } for p in s.parts],
            "overlaps": s.overlaps, "also_in_window": s.also,
            "hints": s.hints, "conflict": s.conflict,
        })
    return {
        "tool": f"peakassign {__version__}",
        "solvent": result.solvent,
        "normalisation": {"raw_integral_per_h": result.unit,
                          "source": result.scale_source},
        "unknown_alert": result.unknown_alert,
        "unknown": [{"ppm": s.ppm, "h": round(s.left_h, 3)}
                    for s in result.category("unknown") + result.hidden],
        "product": species(result.product),
        "reactants": [species(r) for r in result.reactants],
        "declared_impurities": [species(r) for r in result.impurities],
        "peaks": peaks,
        "warnings": result.warnings,
    }


def to_json(result):
    return json.dumps(to_dict(result), ensure_ascii=False, indent=2) + "\n"


CATEGORY_NAME = {
    "product": "产物", "reactant": "原料", "solvent": "溶剂/杂质",
    "unknown": "未知 !!", "ambiguous": "冲突/重叠", "minor": "小峰",
}


def to_text(result):
    """Short summary for the terminal."""
    lines = []
    for state in result.peaks:
        who = "; ".join(p.name if p.role == "solvent"
                        else f"{p.name} {' + '.join(p.labels)}"
                        for p in state.parts)
        lines.append(f"{state.ppm:7.2f}  {state.h:5.2f}H  "
                     f"{CATEGORY_NAME[state.category]:<8}  {who}")
    for report in result.reactants:
        if report.status == "detected":
            flag = "，原料残留较多" if report.alert else ""
            lines.append(f"{report.species.name}: 约为产物的 "
                         f"{report.ratio * 100:.0f} mol%{flag}")
    if result.unknown_alert:
        where = "、".join(f"{s.ppm:.2f} ppm（{s.left_h:.2f}H）"
                         for s in result.category("unknown") + result.hidden)
        lines.append(f"[!] {ALERT}：{where}")
    else:
        lines.append("没有发现积分显著的未知信号。")
    return "\n".join(lines)
