from __future__ import annotations

import csv
import os
import tempfile
from dataclasses import dataclass
from datetime import date
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ANALYSIS_START = date(1993, 10, 1)
ANALYSIS_END = date(2023, 9, 30)
BASE_DIR = Path(__file__).resolve().parent
FLOW_FILE = BASE_DIR / "ID_37 rennsli.csv"
FIGURE_FILE = BASE_DIR / "langaeislina_rennsli.png"


@dataclass
class FlowSeries:
    dates: list[date]
    flow: np.ndarray


@dataclass
class FlowDurationCurve:
    exceedance: np.ndarray
    sorted_flow: np.ndarray
    q5: float
    q50: float
    q95: float

    @property
    def q5_q95_ratio(self) -> float:
        return self.q5 / self.q95

    @property
    def q95_q50_ratio(self) -> float:
        return self.q95 / self.q50


def format_decimal(value: float, decimals: int = 2) -> str:
    return f"{value:.{decimals}f}".replace(".", ",")


def load_flow_series() -> FlowSeries:
    dates: list[date] = []
    flows: list[float] = []

    with FLOW_FILE.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            qobs = row["qobs"].strip()
            if not qobs:
                continue

            current_date = date(int(row["YYYY"]), int(row["MM"]), int(row["DD"]))
            if not (ANALYSIS_START <= current_date <= ANALYSIS_END):
                continue

            dates.append(current_date)
            flows.append(float(qobs))

    if not dates:
        raise ValueError(
            f"Found no valid flow values in {FLOW_FILE} for {ANALYSIS_START.isoformat()} to {ANALYSIS_END.isoformat()}"
        )

    return FlowSeries(dates=dates, flow=np.asarray(flows, dtype=float))


def compute_flow_duration_curve(flow: np.ndarray) -> FlowDurationCurve:
    valid_flow = np.asarray(flow, dtype=float)
    valid_flow = valid_flow[np.isfinite(valid_flow)]
    if len(valid_flow) == 0:
        raise ValueError("No valid flow values available for flow-duration analysis.")

    sorted_flow = np.sort(valid_flow)[::-1]
    exceedance = 100.0 * np.arange(1, len(sorted_flow) + 1, dtype=float) / (len(sorted_flow) + 1)
    q5, q50, q95 = np.interp([5.0, 50.0, 95.0], exceedance, sorted_flow)

    return FlowDurationCurve(
        exceedance=exceedance,
        sorted_flow=sorted_flow,
        q5=float(q5),
        q50=float(q50),
        q95=float(q95),
    )


def plot_flow_duration_curve(curve: FlowDurationCurve) -> None:
    positive = curve.sorted_flow > 0.0

    fig, ax = plt.subplots(figsize=(10, 6.5), constrained_layout=True)
    ax.semilogy(
        curve.exceedance[positive],
        curve.sorted_flow[positive],
        color="#005f73",
        linewidth=2.3,
        label="Langæislína",
    )

    marker_points = [
        (5.0, curve.q5, "Q5", "#ae2012"),
        (50.0, curve.q50, "Q50", "#ee9b00"),
        (95.0, curve.q95, "Q95", "#0a9396"),
    ]
    for exceedance, flow_value, label, color in marker_points:
        ax.axvline(exceedance, color=color, linewidth=1.0, linestyle="--", alpha=0.8)
        ax.scatter(exceedance, flow_value, color=color, s=44, zorder=3)
        ax.annotate(
            f"{label} = {format_decimal(flow_value, 1)} m3/s",
            (exceedance, flow_value),
            xytext=(8, 6),
            textcoords="offset points",
            color=color,
            fontsize=9,
        )

    ax.set_xlim(0, 100)
    ax.set_xlabel("Yfirstigslíkur [% af tíma]")
    ax.set_ylabel("Rennsli [m3/s]")
    ax.set_title(
        f"Langæislína rennslis fyrir Kljáfoss (Hvítá), "
        f"Greiningartímabil: {ANALYSIS_START.isoformat()} til {ANALYSIS_END.isoformat()}"
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper right")

    fig.savefig(FIGURE_FILE, dpi=300)
    plt.close(fig)


def main() -> None:
    series = load_flow_series()
    curve = compute_flow_duration_curve(series.flow)
    plot_flow_duration_curve(curve)
    print(f"Mynd vistuð í {FIGURE_FILE}")


if __name__ == "__main__":
    main()
