from __future__ import annotations

import csv
import math
import os
import tempfile
from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CATCHMENT_ID = "37"
ALPHA = 0.925
LADSON_PASSES = 3
LADSON_REFLECTION = 30
SENSITIVITY_PASSES = 5
ANALYSIS_START = date(1993, 10, 1)
ANALYSIS_END = date(2023, 9, 30)
REC_MAX_DAILY_PRECIP = 1.0
REC_MAX_3DAY_PRECIP = 2.0
REC_MIN_LENGTH = 7
REC_MIN_R2 = 0.90
REC_MAX_QC_FLAG = 80

BASE_DIR = Path(__file__).resolve().parent
FLOW_FILE = BASE_DIR / "ID_37 rennsli.csv"
WEATHER_FILE = BASE_DIR / "ID_37_rett.csv"
WEATHER_FALLBACK_FILE = BASE_DIR.parent / "ID_37.csv"
CATCHMENT_FILE = BASE_DIR.parent / "Catchment_attributes.csv"
GAUGE_FILE = BASE_DIR.parent / "lamah_ice" / "D_gauges" / "1_attributes" / "Gauge_attributes.csv"
FIGURE_FILE = BASE_DIR / "grunnvatnsframlag_greining.png"
SUMMARY_FILE = BASE_DIR / "grunnvatnsframlag_nidurstodur.txt"


@dataclass
class FlowSeries:
    # Heldur utan um samkeyrða tíma-, rennslis- og veðurröð fyrir stöðina.
    dates: list[date]
    flow: np.ndarray
    quality: np.ndarray
    precipitation: np.ndarray
    temperature: np.ndarray


@dataclass
class RecessionSegment:
    start_index: int
    end_index: int
    slope: float
    intercept: float
    r_squared: float

    @property
    def length(self) -> int:
        return self.end_index - self.start_index + 1


@dataclass
class RecessionResult:
    segments: list[RecessionSegment]
    master_offsets: np.ndarray
    master_curve: np.ndarray
    master_counts: np.ndarray
    slope: float
    intercept: float

    @property
    def recession_ratio(self) -> float:
        return math.exp(self.slope)

    @property
    def recession_constant(self) -> float:
        return -self.slope

    @property
    def recession_timescale_days(self) -> float:
        return -1.0 / self.slope


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        return float(text.replace(",", "."))
    except ValueError:
        # Sum afrit af veðurgögnunum hafa verið vistuð í töflureikni og fengið dagsetningasnið.
        return None


def load_flow_series() -> FlowSeries:
    dates: list[date] = []
    flows: list[float] = []
    quality: list[int] = []

    # Les daglegt rennsli og gæðaflögg stöðvar 37.
    with FLOW_FILE.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            dates.append(date(int(row["YYYY"]), int(row["MM"]), int(row["DD"])))
            flows.append(float(row["qobs"]))
            quality.append(int(float(row["qc_flag"])))

    weather_by_date: dict[date, tuple[float | None, float | None]] = {}
    # Les úrkomu og hita svo hægt sé að velja þurr recession-skeið.
    weather_path = WEATHER_FALLBACK_FILE if WEATHER_FALLBACK_FILE.exists() else WEATHER_FILE
    with weather_path.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            current_date = date(int(row["YYYY"]), int(row["MM"]), int(row["DD"]))
            weather_by_date[current_date] = (
                parse_float(row["prec"]),
                parse_float(row["2m_temp_mean"]),
            )

    precipitation = []
    temperature = []
    for current_date in dates:
        prec, temp = weather_by_date.get(current_date, (None, None))
        precipitation.append(np.nan if prec is None else prec)
        temperature.append(np.nan if temp is None else temp)

    analysis = period_mask(dates, ANALYSIS_START, ANALYSIS_END)

    return FlowSeries(
        dates=[current_date for current_date, keep in zip(dates, analysis) if keep],
        flow=np.asarray(flows, dtype=float)[analysis],
        quality=np.asarray(quality, dtype=int)[analysis],
        precipitation=np.asarray(precipitation, dtype=float)[analysis],
        temperature=np.asarray(temperature, dtype=float)[analysis],
    )


def load_csv_row(path: Path, row_id: str) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            if row["id"] == row_id:
                return row
    raise ValueError(f"Found no row with id={row_id} in {path}")


def lyne_hollick_baseflow(
    flow: np.ndarray,
    alpha: float = ALPHA,
    passes: int = LADSON_PASSES,
    reflection: int = LADSON_REFLECTION,
) -> np.ndarray:
    if len(flow) < reflection + 2:
        raise ValueError("Time series is too short for the requested reflection length.")

    # Speglum raðirnar í báða enda til að draga úr "warm-up" áhrifum síunnar.
    padded = np.concatenate(
        [
            flow[reflection:0:-1],
            flow,
            flow[-2 : -reflection - 2 : -1],
        ]
    )
    n_values = len(padded)
    starts = ([0, n_values - 1] * passes)[:passes]
    ends = ([n_values - 1, 0] * passes)[:passes]
    increments = ([1, -1] * passes)[:passes]

    quickflow = np.zeros(n_values, dtype=float)
    quickflow[0] = padded[0]
    baseflow = padded.copy()

    # Þetta er stöðluð Ladson/Lyne-Hollick síun þar sem quickflow er fjarlægt í hverri yfirferð.
    for pass_index in range(passes):
        for i in range(starts[pass_index] + increments[pass_index], ends[pass_index] + increments[pass_index], increments[pass_index]):
            quickflow[i] = (
                alpha * quickflow[i - increments[pass_index]]
                + ((1.0 + alpha) / 2.0) * (baseflow[i] - baseflow[i - increments[pass_index]])
            )
        positive = quickflow > 0
        baseflow[positive] = baseflow[positive] - quickflow[positive]
        quickflow[ends[pass_index]] = baseflow[ends[pass_index]]

    return baseflow[reflection:-reflection]


def baseflow_index(flow: np.ndarray, baseflow: np.ndarray) -> float:
    return float(np.sum(baseflow) / np.sum(flow))


def period_mask(dates: list[date], start: date, end: date) -> np.ndarray:
    return np.asarray([start <= current_date <= end for current_date in dates], dtype=bool)


def monthly_climatology(dates: list[date], values: np.ndarray) -> np.ndarray:
    totals = np.zeros(12, dtype=float)
    counts = np.zeros(12, dtype=int)
    for current_date, value in zip(dates, values):
        totals[current_date.month - 1] += value
        counts[current_date.month - 1] += 1
    return totals / counts


def daily_water_year_climatology(dates: list[date], values: np.ndarray) -> tuple[list[date], np.ndarray]:
    # Notum eitt viðmiðunar-vatnsár til að byggja daglegt meðaltalsár frá október til september.
    reference_dates: list[date] = []
    current_date = date(2001, 10, 1)
    end_date = date(2002, 9, 30)
    while current_date <= end_date:
        reference_dates.append(current_date)
        current_date = date.fromordinal(current_date.toordinal() + 1)

    grouped_values: dict[tuple[int, int], list[float]] = {(d.month, d.day): [] for d in reference_dates}
    for current_date, value in zip(dates, values):
        if current_date.month == 2 and current_date.day == 29:
            continue
        grouped_values[(current_date.month, current_date.day)].append(float(value))

    climatology = np.asarray(
        [np.mean(grouped_values[(d.month, d.day)]) for d in reference_dates],
        dtype=float,
    )
    return reference_dates, climatology


def fit_log_linear(x_values: np.ndarray, y_values: np.ndarray) -> tuple[float, float, float]:
    slope, intercept = np.polyfit(x_values, np.log(y_values), deg=1)
    fitted = slope * x_values + intercept
    residuals = np.log(y_values) - fitted
    total = np.log(y_values) - np.mean(np.log(y_values))
    if np.allclose(total, 0.0):
        r_squared = 1.0
    else:
        r_squared = 1.0 - float(np.sum(residuals**2) / np.sum(total**2))
    return float(slope), float(intercept), r_squared


def extract_recession_segments(flow_series: FlowSeries) -> RecessionResult:
    candidate_segments: list[tuple[int, int]] = []
    start_index: int | None = None

    for i in range(1, len(flow_series.flow)):
        recent_precip = flow_series.precipitation[max(0, i - 2) : i + 1]
        # Recession-skeið þurfa að vera þurr, hlý og með lækkandi rennsli.
        is_dry = (
            np.isfinite(flow_series.precipitation[i])
            and flow_series.precipitation[i] <= REC_MAX_DAILY_PRECIP
            and float(np.nansum(recent_precip)) <= REC_MAX_3DAY_PRECIP
        )
        is_warm = np.isfinite(flow_series.temperature[i]) and flow_series.temperature[i] > 0.0
        has_good_quality = (
            flow_series.quality[i] <= REC_MAX_QC_FLAG
            and flow_series.quality[i - 1] <= REC_MAX_QC_FLAG
        )
        is_receding = flow_series.flow[i] <= flow_series.flow[i - 1]

        if is_dry and is_warm and has_good_quality and is_receding:
            if start_index is None:
                start_index = i - 1
        else:
            if start_index is not None and i - start_index >= REC_MIN_LENGTH:
                candidate_segments.append((start_index, i - 1))
            start_index = None

    if start_index is not None and len(flow_series.flow) - start_index >= REC_MIN_LENGTH:
        candidate_segments.append((start_index, len(flow_series.flow) - 1))

    accepted_segments: list[RecessionSegment] = []
    for start, end in candidate_segments:
        x_values = np.arange(end - start + 1, dtype=float)
        y_values = flow_series.flow[start : end + 1]
        slope, intercept, r_squared = fit_log_linear(x_values, y_values)
        # Við höldum aðeins skeiðum sem líkjast vel veldisvísisrýrnun.
        if slope < 0.0 and r_squared >= REC_MIN_R2:
            accepted_segments.append(
                RecessionSegment(
                    start_index=start,
                    end_index=end,
                    slope=slope,
                    intercept=intercept,
                    r_squared=r_squared,
                )
            )

    if not accepted_segments:
        raise ValueError("No recession segments satisfied the selected criteria.")

    max_length = max(segment.length for segment in accepted_segments)
    master_curve = []
    master_counts = []
    # Byggjum miðgildis rýrnunarferil úr miðgildi staðlaðra skeiða Q/Q0.
    for offset in range(max_length):
        values_at_offset = [
            flow_series.flow[segment.start_index + offset] / flow_series.flow[segment.start_index]
            for segment in accepted_segments
            if offset < segment.length
        ]
        master_curve.append(float(np.median(values_at_offset)))
        master_counts.append(len(values_at_offset))

    master_curve_array = np.asarray(master_curve, dtype=float)
    master_counts_array = np.asarray(master_counts, dtype=int)
    valid = master_counts_array >= 5
    master_offsets = np.arange(max_length, dtype=float)[valid]
    master_log = np.log(master_curve_array[valid])

    # Hallatala í log-rúmi gefur a, og k fæst síðan sem exp(-a).
    slope = float(np.dot(master_offsets, master_log) / np.dot(master_offsets, master_offsets))
    intercept = 0.0

    return RecessionResult(
        segments=accepted_segments,
        master_offsets=master_offsets,
        master_curve=master_curve_array[valid],
        master_counts=master_counts_array[valid],
        slope=slope,
        intercept=intercept,
    )


def write_summary(
    flow_series: FlowSeries,
    gauge: dict[str, str],
    catchment: dict[str, str],
    period_bfi: float,
    sensitivity_bfi: float,
    recession: RecessionResult,
) -> str:
    quality_counts = Counter(flow_series.quality.tolist())
    total_count = len(flow_series.flow)
    quality_text = ", ".join(
        f"qc={flag}: {count} dagar ({100.0 * count / total_count:.1f}%)"
        for flag, count in sorted(quality_counts.items())
    )

    # Útbúum texta sem hægt er að nota beint í verkefnaskil.
    summary = "\n".join(
        [
            f"Mat á grunnvatnsframlagi fyrir stöð {CATCHMENT_ID} ({gauge['name']} í {gauge['river']})",
            "",
            "Gögn:",
            f"- Rennsli: {flow_series.dates[0].isoformat()} til {flow_series.dates[-1].isoformat()} ({total_count} dagar)",
            f"- Vatnasviðsflatarmál: {float(catchment['area_calc']):.1f} km2",
            "",
            "Aðferð:",
            f"- Baseflow separation með Lyne-Hollick/Ladson (alpha={ALPHA}, passes={LADSON_PASSES}, reflection={LADSON_REFLECTION})",
            "- Rýrnunargreining a þurrum og hallandi rennslisskeiðum (P <= 1 mm/d, 3 daga summa <= 2 mm, T > 0 C, qc <= 80, mín. 7 dagar, R2 >= 0.90)",
            "",
            "Helstu niðurstöður:",
            f"- BFI fyrir tímabilið {ANALYSIS_START.isoformat()} til {ANALYSIS_END.isoformat()} med 3 ferlum: {period_bfi:.3f}",
            f"- Næmnipróf med 5 ferlum yfir sama tímabil: {sensitivity_bfi:.3f}",
            f"- Daglegt rýrnunarhlutfall k (Q(t+1)/Q(t) ur miðgildisferli): {recession.recession_ratio:.3f} [-]",
            f"- Rýrnunarstuðull a i Q(t)=Q0*exp(-a*t): {recession.recession_constant:.4f} d^-1",
            f"- Einkennandi tæmingartími : {recession.recession_timescale_days:.0f} dagar",
            f"- Fjöldi samþykktra rýrnunar-skeiða: {len(recession.segments)}",
        ]
    )

    SUMMARY_FILE.write_text(summary + "\n", encoding="utf-8")
    return summary


def plot_results(
    flow_series: FlowSeries,
    baseflow: np.ndarray,
    recession: RecessionResult,
    gauge: dict[str, str],
    period_bfi: float,
) -> None:
    month_labels = ["okt", "nov", "des", "jan", "feb", "mar", "apr", "mai", "jun", "jul", "agu", "sep"]
    hydrological_order = np.array([10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]) - 1
    # Efsta grafið sýnir nú meðaldagsferil yfir eitt vatnsár.
    plot_dates, plot_flow = daily_water_year_climatology(flow_series.dates, flow_series.flow)
    _, plot_baseflow = daily_water_year_climatology(flow_series.dates, baseflow)
    monthly_flow = monthly_climatology(flow_series.dates, flow_series.flow)[hydrological_order]
    monthly_baseflow = monthly_climatology(flow_series.dates, baseflow)[hydrological_order]
    month_ticks = [
        date(2001, 10, 1),
        date(2001, 11, 1),
        date(2001, 12, 1),
        date(2002, 1, 1),
        date(2002, 2, 1),
        date(2002, 3, 1),
        date(2002, 4, 1),
        date(2002, 5, 1),
        date(2002, 6, 1),
        date(2002, 7, 1),
        date(2002, 8, 1),
        date(2002, 9, 1),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(13, 14), constrained_layout=True)

    axes[0].fill_between(plot_dates, 0.0, plot_baseflow, color="#8ecae6", alpha=0.9, label="Grunnrennsli (Baseflow)")
    axes[0].fill_between(plot_dates, plot_baseflow, plot_flow, color="#ffb703", alpha=0.65, label="Yfirborðsrennsli (Quickflow)")
    axes[0].plot(plot_dates, plot_flow, color="#023047", linewidth=1.2, label="Heildarrennsli")
    axes[0].set_xticks(month_ticks)
    axes[0].set_xticklabels(month_labels)
    axes[0].set_ylabel("Rennsli [m3/s]")
    axes[0].set_title(
        f"Meðalrennslisferill eins vatnsárs fyrir {gauge['name']} ({gauge['river']})\n"
        f"BFI tímabilsins {ANALYSIS_START.isoformat()} til {ANALYSIS_END.isoformat()} = {period_bfi:.3f}"
    )
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(loc="upper right")

    month_numbers = np.arange(1, 13)
    axes[1].plot(month_numbers, monthly_flow, color="#023047", linewidth=2.0, marker="o", label="Heildarrennsli")
    axes[1].plot(month_numbers, monthly_baseflow, color="#219ebc", linewidth=2.0, marker="o", label="Grunnrennsli (Baseflow)")
    axes[1].set_xticks(month_numbers)
    axes[1].set_xticklabels(month_labels)
    axes[1].set_ylabel("Meðalrennsli [m3/s]")
    axes[1].set_title("Mánaðarlegt meðaltal rennslis og grunnrennslis eftir vatnsári")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(loc="upper right")

    # Gráu línurnar sýna einstök rýrnunar-skeið, blá línan miðgildi þeirra.
    for segment in recession.segments:
        offsets = np.arange(segment.length)
        normalized = flow_series.flow[segment.start_index : segment.end_index + 1] / flow_series.flow[segment.start_index]
        axes[2].semilogy(offsets, normalized, color="#adb5bd", linewidth=0.9, alpha=0.85)

    fitted_curve = np.exp(recession.intercept + recession.slope * recession.master_offsets)
    axes[2].semilogy(
        recession.master_offsets,
        recession.master_curve,
        color="#005f73",
        linewidth=2.5,
        marker="o",
        label="Miðgildis rýrnunarferill",
    )
    axes[2].semilogy(
        recession.master_offsets,
        fitted_curve,
        color="#ae2012",
        linewidth=2.0,
        linestyle="--",
        label=f"Aðhvarfslína: Q/Q0 = exp({recession.slope:.3f} t)",
    )
    axes[2].set_xlabel("Dagar frá byrjun rýrnunar")
    axes[2].set_ylabel("Q / Q0 [-]")
    axes[2].set_title(
        f"Rýrnunargreining ({len(recession.segments)} ásættanleg skeið, "
        f"k={recession.recession_ratio:.3f}, a={recession.recession_constant:.4f} d^-1)"
    )
    axes[2].grid(True, which="both", alpha=0.25)
    axes[2].legend(loc="upper right")

    fig.savefig(FIGURE_FILE, dpi=300)
    plt.close(fig)


def main() -> None:
    flow_series = load_flow_series()
    baseflow = lyne_hollick_baseflow(flow_series.flow)
    period_bfi = baseflow_index(flow_series.flow, baseflow)

    # Næmnipróf á sama greiningartímabili með fleiri filter passes.
    sensitivity_baseflow = lyne_hollick_baseflow(flow_series.flow, passes=SENSITIVITY_PASSES)
    sensitivity_bfi = baseflow_index(flow_series.flow, sensitivity_baseflow)

    recession = extract_recession_segments(flow_series)
    gauge = load_csv_row(GAUGE_FILE, CATCHMENT_ID)
    catchment = load_csv_row(CATCHMENT_FILE, CATCHMENT_ID)

    summary = write_summary(
        flow_series=flow_series,
        gauge=gauge,
        catchment=catchment,
        period_bfi=period_bfi,
        sensitivity_bfi=sensitivity_bfi,
        recession=recession,
    )
    plot_results(
        flow_series=flow_series,
        baseflow=baseflow,
        recession=recession,
        gauge=gauge,
        period_bfi=period_bfi,
    )

    print(summary)
    print(f"\nMynd vistud i {FIGURE_FILE}")
    print(f"Nidurstodur vistaðar i {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
