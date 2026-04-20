import csv
import os
import tempfile
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


START_DATE = date(1993, 10, 1)
END_DATE = date(2023, 9, 30)
BASE_DIR = Path(__file__).resolve().parent
FLOW_FILE = BASE_DIR / "ID_37 rennsli.csv"
WEATHER_FILE = BASE_DIR / "ID_37_rett.csv"
OUTPUT_FILE = BASE_DIR / "medaltalsar_arstidarsveifla.png"
TEMP_COLUMN = "2m_temp_carra"
PREC_COLUMN = "prec_carra"


def parse_float(value: str) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    return float(text)


def in_period(current_date: date) -> bool:
    return START_DATE <= current_date <= END_DATE


def hydrological_plot_date(month: int) -> datetime:
    year = 1999 if month >= 10 else 2000
    return datetime(year, month, 1)


def format_date_range(dates: list[date]) -> str:
    if not dates:
        return "engin gildi fundust"
    return f"{min(dates).isoformat()} til {max(dates).isoformat()}"


def read_flow_means() -> list[tuple[datetime, float]]:
    grouped_values: dict[int, list[float]] = defaultdict(list)

    with FLOW_FILE.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            current_date = date(int(row["YYYY"]), int(row["MM"]), int(row["DD"]))
            if not in_period(current_date):
                continue

            flow = parse_float(row["qobs"])
            if flow is not None:
                grouped_values[current_date.month].append(flow)

    means = []
    for month, values in grouped_values.items():
        means.append((hydrological_plot_date(month), sum(values) / len(values)))
    means.sort(key=lambda item: item[0])
    return means


def read_weather_means() -> tuple[list[tuple[datetime, float]], list[tuple[datetime, float]]]:
    temp_values: dict[int, list[float]] = defaultdict(list)
    prec_values: dict[int, list[float]] = defaultdict(list)
    temp_dates_used: list[date] = []
    prec_dates_used: list[date] = []

    with WEATHER_FILE.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            current_date = date(int(row["YYYY"]), int(row["MM"]), int(row["DD"]))
            if not in_period(current_date):
                continue

            temp = parse_float(row[TEMP_COLUMN])
            prec = parse_float(row[PREC_COLUMN])

            if temp is not None:
                temp_values[current_date.month].append(temp)
                temp_dates_used.append(current_date)
            if prec is not None:
                prec_values[current_date.month].append(prec)
                prec_dates_used.append(current_date)

    temp_means = []
    for month, values in temp_values.items():
        temp_means.append((hydrological_plot_date(month), sum(values) / len(values)))
    temp_means.sort(key=lambda item: item[0])

    prec_means = []
    for month, values in prec_values.items():
        prec_means.append((hydrological_plot_date(month), sum(values) / len(values)))
    prec_means.sort(key=lambda item: item[0])

    print(f"Notad hitastig ur {TEMP_COLUMN}: {format_date_range(temp_dates_used)}")
    print(f"Notud urkoma ur {PREC_COLUMN}: {format_date_range(prec_dates_used)}")

    return temp_means, prec_means


def month_ticks() -> tuple[list[datetime], list[str]]:
    ticks = [
        datetime(1999, 10, 1),
        datetime(1999, 11, 1),
        datetime(1999, 12, 1),
        datetime(2000, 1, 1),
        datetime(2000, 2, 1),
        datetime(2000, 3, 1),
        datetime(2000, 4, 1),
        datetime(2000, 5, 1),
        datetime(2000, 6, 1),
        datetime(2000, 7, 1),
        datetime(2000, 8, 1),
        datetime(2000, 9, 1),
    ]
    labels = ["okt", "nov", "des", "jan", "feb", "mar", "apr", "maí", "jún", "júl", "ág", "sep"]
    return ticks, labels


def plot_mean_annual_cycle(
    flow_means: list[tuple[datetime, float]],
    temp_means: list[tuple[datetime, float]],
    prec_means: list[tuple[datetime, float]],
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True, constrained_layout=True)

    flow_dates = [item[0] for item in flow_means]
    flow_values = [item[1] for item in flow_means]
    prec_dates = [item[0] for item in prec_means]
    prec_values = [item[1] for item in prec_means]
    temp_dates = [item[0] for item in temp_means]
    temp_values = [item[1] for item in temp_means]

    axes[0].plot(
        flow_dates,
        flow_values,
        color="#0843C3",
        linewidth=2,
        marker="o",
        markersize=5,
    )
    axes[0].set_ylabel("Rennsli [m3/s]")
    axes[0].set_title("Meðaltalsár tímabilsins 1.10.1993 til 30.9.2023")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(
        prec_dates,
        prec_values,
        color="#2b8a03",
        linewidth=2,
        marker="o",
        markersize=5,
    )
    axes[1].set_ylabel("Úrkoma CARRA [mm/d]")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(
        temp_dates,
        temp_values,
        color="#bb3e03",
        linewidth=2,
        marker="o",
        markersize=5,
    )
    axes[2].axhline(0, color="black", linewidth=0.9, alpha=0.6)
    axes[2].set_ylabel("Hitastig CARRA [°C]")
    axes[2].set_xlabel("Vatnsár (október til september)")
    axes[2].grid(True, alpha=0.3)

    ticks, labels = month_ticks()
    axes[2].set_xticks(ticks)
    axes[2].set_xticklabels(labels)

    fig.savefig(OUTPUT_FILE, dpi=300)
    plt.close(fig)


def main() -> None:
    flow_means = read_flow_means()
    temp_means, prec_means = read_weather_means()
    plot_mean_annual_cycle(flow_means, temp_means, prec_means)
    print(f"Mynd vistud i {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
