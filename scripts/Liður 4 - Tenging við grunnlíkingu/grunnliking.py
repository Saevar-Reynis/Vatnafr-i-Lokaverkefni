from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import date
from pathlib import Path


ANALYSIS_START = date(1993, 10, 1)
ANALYSIS_END = date(2023, 9, 30)
CATCHMENT_ID = "37"

BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent

FLOW_FILE = BASE_DIR / "ID_37 rennsli.csv"
WEATHER_FILE = BASE_DIR / "ID_37_rett.csv"
CATCHMENT_FILE = PROJECT_ROOT / "Catchment_attributes.csv"

TABLE_MD_FILE = BASE_DIR / "grunnliking_tafla.md"
TABLE_CSV_FILE = BASE_DIR / "grunnliking_tafla.csv"


@dataclass
class WaterYearTotals:
    days: int = 0
    p_mm: float = 0.0
    q_mm: float = 0.0
    et_mm: float = 0.0
    q_m3s_sum: float = 0.0


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    return float(text)


def load_csv_row(path: Path, row_id: str) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            if row["id"] == row_id:
                return row
    raise ValueError(f"Found no row with id={row_id} in {path}")


def water_year(current_date: date) -> int:
    return current_date.year + 1 if current_date.month >= 10 else current_date.year


def runoff_depth_mm_per_day(flow_m3s: float, area_km2: float) -> float:
    return flow_m3s * 86.4 / area_km2


def compute_statistics() -> tuple[float, float, float, float]:
    catchment = load_csv_row(CATCHMENT_FILE, CATCHMENT_ID)
    area_km2 = float(catchment["area_calc"])

    weather_by_date: dict[date, tuple[float, float]] = {}
    with WEATHER_FILE.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            current_date = date(
                int(float(row["YYYY"])),
                int(float(row["MM"])),
                int(float(row["DD"])),
            )
            if not (ANALYSIS_START <= current_date <= ANALYSIS_END):
                continue

            prec_carra = parse_float(row["prec_carra"])
            total_et_carra = parse_float(row["total_et_carra"])
            if prec_carra is not None and total_et_carra is not None:
                weather_by_date[current_date] = (prec_carra, total_et_carra)

    totals_by_year: dict[int, WaterYearTotals] = {}
    with FLOW_FILE.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            current_date = date(int(row["YYYY"]), int(row["MM"]), int(row["DD"]))
            if not (ANALYSIS_START <= current_date <= ANALYSIS_END):
                continue
            if current_date not in weather_by_date:
                continue

            wy = water_year(current_date)
            if wy not in totals_by_year:
                totals_by_year[wy] = WaterYearTotals()

            q_m3s = float(row["qobs"])
            prec_carra, total_et_carra = weather_by_date[current_date]
            totals = totals_by_year[wy]
            totals.days += 1
            totals.p_mm += prec_carra
            totals.q_mm += runoff_depth_mm_per_day(q_m3s, area_km2)
            totals.et_mm += total_et_carra
            totals.q_m3s_sum += q_m3s

    full_years = [totals for _, totals in sorted(totals_by_year.items()) if totals.days >= 365]
    if not full_years:
        raise ValueError("No full water years found in the selected period.")

    mean_p_mm = sum(item.p_mm for item in full_years) / len(full_years)
    mean_q_mm = sum(item.q_mm for item in full_years) / len(full_years)
    mean_et_mm = sum(item.et_mm for item in full_years) / len(full_years)
    mean_q_m3s = sum(item.q_m3s_sum / item.days for item in full_years) / len(full_years)

    return mean_p_mm, mean_q_mm, mean_et_mm, mean_q_m3s


def build_table_rows(mean_p_mm: float, mean_q_mm: float, mean_et_mm: float, mean_q_m3s: float) -> list[list[str]]:
    return [
        ["P úrkoma", "CARRA / vatnasviðsmeðaltal", f"{mean_p_mm:.0f} mm/ári"],
        ["Q afrennsli", "mælt rennsli", f"{mean_q_mm:.0f} mm/ári eða {mean_q_m3s:.1f} m³/s"],
        ["ET", "CARRA / total_et_carra", f"{mean_et_mm:.0f} mm/ári"],
        ["ΔS", "ekki mælt beint", "—"],
        ["G", "ekki mælt", "—"],
    ]


def write_table_files(rows: list[list[str]]) -> None:
    headers = ["Þáttur", "Hvernig metinn", "Stærð yfir tímabilið"]

    with TABLE_CSV_FILE.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(rows)

    markdown_lines = [
        "| Þáttur | Hvernig metinn | Stærð yfir tímabilið |",
        "|---|---|---|",
    ]
    for row in rows:
        markdown_lines.append(f"| {row[0]} | {row[1]} | {row[2]} |")
    TABLE_MD_FILE.write_text("\n".join(markdown_lines) + "\n", encoding="utf-8")


def main() -> None:
    mean_p_mm, mean_q_mm, mean_et_mm, mean_q_m3s = compute_statistics()
    table_rows = build_table_rows(mean_p_mm, mean_q_mm, mean_et_mm, mean_q_m3s)
    write_table_files(table_rows)

    print(f"Tafla vistuð í {TABLE_MD_FILE}")
    print(f"Tafla vistuð í {TABLE_CSV_FILE}")


if __name__ == "__main__":
    main()
