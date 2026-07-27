from __future__ import annotations

import csv
import json
import re
from pathlib import Path


VIEWER_DIR = Path(__file__).resolve().parent
RESULT_ROOT = Path(
    r"D:\Dropbox_new\Dropbox\2_저널페이퍼\01_cscm_blast\05_dawon_2d_beam"
    r"\01_analysis_2dbeam_blast\251222_dawon_m0.1_to_m0.8_copy2"
)
OUTPUT_DIR = VIEWER_DIR / "data" / "structure_damage"

COLS = 220
ROWS = 8
DY_CM = 0.5
DZ_CM = 0.5
Y0_CM = -55.0
Z0_CM = -4.0


def charge_tag_to_folder(tag: str) -> str:
    return "m=" + tag.removeprefix("m").replace("p", ".")


def case_file_tag(charge_tag: str, x_cm: float, y_cm: float) -> str:
    return f"{charge_tag}_x{int(round(x_cm)):03d}_y{int(round(y_cm)):03d}.json"


def parse_failed_elements(case_dir: Path) -> set[int]:
    failed: set[int] = set()
    log_path = case_dir / "messag"
    if not log_path.exists():
        return failed

    pattern = re.compile(r"solid element\s+(\d+)\s+failed")
    with log_path.open(errors="ignore") as file:
        for line in file:
            match = pattern.search(line)
            if match:
                failed.add(int(match.group(1)))
    return failed


def parse_damage_history(case_dir: Path) -> tuple[list[float], int]:
    csv_path = case_dir / "damage_history_all_elements.csv"
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)

    with csv_path.open(newline="", errors="ignore") as file:
        _title = file.readline()
        header = file.readline().strip().split(",")
        element_ids = [int(value) for value in header[1:] if value.strip()]
        maxima = [0.0] * len(element_ids)

        for line in file:
            if not line.strip():
                continue
            tokens = line.split(",")
            for index in range(len(element_ids)):
                token_index = index + 1
                if token_index >= len(tokens):
                    continue
                token = tokens[token_index].strip()
                if not token:
                    continue
                try:
                    value = float(token)
                except ValueError:
                    continue
                if value > maxima[index]:
                    maxima[index] = value

    failed = parse_failed_elements(case_dir)
    failed_count = 0
    for index, element_id in enumerate(element_ids):
        if element_id in failed:
            maxima[index] = 1.0
            failed_count += 1

    return maxima, failed_count


def write_geometry() -> None:
    geometry = {
        "cols": COLS,
        "rows": ROWS,
        "y_centers_cm": [round(Y0_CM + (index + 0.5) * DY_CM, 3) for index in range(COLS)],
        "z_centers_cm": [round(Z0_CM + (index + 0.5) * DZ_CM, 3) for index in range(ROWS)],
        "value_order": "for each y column, z rows progress from bottom to top",
        "averaging": "area-weighted over selected structural concrete elements; eroded elements are assigned d=1.0",
    }
    (OUTPUT_DIR / "geometry.json").write_text(json.dumps(geometry, separators=(",", ":")), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_geometry()

    manifest = {}
    missing = []
    for heatmap_csv in sorted((VIEWER_DIR / "data").glob("m*.csv")):
        charge_tag = heatmap_csv.stem
        charge_folder = charge_tag_to_folder(charge_tag)
        manifest[charge_tag] = {}

        with heatmap_csv.open(newline="", encoding="utf-8-sig") as file:
            for row in csv.DictReader(file):
                x_cm = float(row["x_cm"])
                y_cm = float(row["y_cm"])
                case_dir = RESULT_ROOT / charge_folder / "left" / "v" / f"y{int(round(x_cm)):03d}_z{int(round(y_cm)):03d}"
                output_name = case_file_tag(charge_tag, x_cm, y_cm)
                output_path = OUTPUT_DIR / output_name

                try:
                    values, failed_count = parse_damage_history(case_dir)
                except FileNotFoundError:
                    missing.append(str(case_dir))
                    continue

                if len(values) != COLS * ROWS:
                    raise ValueError(f"{case_dir} has {len(values)} elements, expected {COLS * ROWS}")

                payload = {
                    "charge": charge_tag,
                    "x_cm": x_cm,
                    "y_cm": y_cm,
                    "failed_elements": failed_count,
                    "values": [round(value, 6) for value in values],
                }
                output_path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
                manifest[charge_tag][f"{int(round(x_cm)):03d}_{int(round(y_cm)):03d}"] = output_name

    (OUTPUT_DIR / "manifest.json").write_text(json.dumps(manifest, separators=(",", ":")), encoding="utf-8")
    if missing:
        (OUTPUT_DIR / "missing_cases.txt").write_text("\n".join(missing), encoding="utf-8")
        print(f"Missing cases: {len(missing)}")
    print(f"Wrote structural damage files to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
