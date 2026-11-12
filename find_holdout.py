from __future__ import annotations

import json
import random
import re
import subprocess
from pathlib import Path

import pyarrow.parquet as pq

PARQUET = Path(r"C:\Users\YINGJIE SHAN\Desktop\NMR project\NMRexp_10to24_1_1004.parquet")
REPORTS_DIR = Path("reports")
SAMPLE_SIZE = 40
RANDOM_SEED = 42


def load_smiles_column(parquet_path: Path) -> list[str]:
    schema = pq.read_schema(parquet_path)
    candidates = [n for n in schema.names if "smiles" in n.lower()]
    if not candidates:
        raise SystemExit(f"在 {parquet_path} 里找不到 SMILES 列，实际列名：{schema.names}")
    col = candidates[0]
    table = pq.read_table(parquet_path, columns=[col])
    values = table.column(col).to_pylist()
    seen: set[str] = set()
    out: list[str] = []
    for v in values:
        if not v or not isinstance(v, str):
            continue
        v = v.strip()
        if not v or v in seen:
            continue
        if len(v) < 4 or " " in v:
            continue
        seen.add(v)
        out.append(v)
    return out


def run_check(smiles: str) -> dict | None:
    proc = subprocess.run(
        ["nmrcheck", "check", smiles],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if proc.returncode != 0:
        return None
    m = re.search(r"json:\s*(\S+)", proc.stdout)
    if not m:
        return None
    report_path = Path(m.group(1))
    if not report_path.exists():
        alt = REPORTS_DIR / report_path.name
        if alt.exists():
            report_path = alt
        else:
            return None
    try:
        return json.loads(report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def has_independent_warning(report: dict) -> bool:
    for entry in report.get("degradation_log", []) or []:
        code = entry.get("code") if isinstance(entry, dict) else None
        if code == "prediction_not_independent":
            return True
    text = json.dumps(report, ensure_ascii=False).lower()
    return "prediction_not_independent" in text or "not independent" in text


def record_count(report: dict) -> int:
    text = json.dumps(report, ensure_ascii=False)
    return text.count('"experimental"')


def main() -> None:
    if not PARQUET.exists():
        raise SystemExit(f"找不到 parquet：{PARQUET}")

    print(f"读取 {PARQUET.name} 的 SMILES 列 ...")
    all_smiles = load_smiles_column(PARQUET)
    print(f"共 {len(all_smiles)} 个去重后的 SMILES")

    rng = random.Random(RANDOM_SEED)
    sample = rng.sample(all_smiles, min(SAMPLE_SIZE, len(all_smiles)))

    tried = 0
    for i, smiles in enumerate(sample, 1):
        tried += 1
        print(f"[{i}/{len(sample)}] 尝试 {smiles}", end=" ... ", flush=True)
        report = run_check(smiles)
        if report is None:
            print("check 失败，跳过")
            continue

        if has_independent_warning(report):
            print("在训练集里，跳过")
            continue

        n = record_count(report)
        verdict = report.get("dedup_verdict") or report.get("verdict") or "?"
        conf = report.get("confidence", "?")
        print(f"\n找到训练集外化合物：{smiles}")
        print(f"  verdict: {verdict}, confidence: {conf}, experimental 记录约 {n} 条")
        print(f"  请到 reports/ 文件夹下找对应的 .md 报告，确认 Degradation log 里没有 prediction_not_independent")
        return

    print(f"\n尝试了 {tried} 个化合物，全部都在预测器训练集里。")
    print("建议：把 SAMPLE_SIZE 调大（比如 200），或换 RANDOM_SEED 重跑。")


if __name__ == "__main__":
    main()