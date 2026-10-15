"""Reporting-only projection of prideQC's opt-in TMT mzQC v1 evidence.

No reporter masses, hits, denominators or percentages are inferred here. All
numeric values are copied from the original mzQC record; stored fractions are
multiplied by 100 solely for presentation. Catalogs and MS levels stay separate.
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

CENTROID_ACCESSION = "QCPRIDE:D9D5612FDC49DA507254"
PROFILE_ACCESSION = "QCPRIDE:94464A800A02DA828D7E"

_TMT_SOURCES = {
    CENTROID_ACCESSION: ("centroid", "EXPLORATORY_NOT_FOR_AUTOMATIC_SDRF_WRITEBACK"),
    PROFILE_ACCESSION: (
        "profile",
        "EXPERIMENTAL_PROFILE_PEAK_PICKING_NOT_FOR_SDRF_WRITEBACK",
    ),
}

# Every numeric column added to the HTML TMT sections has this originating
# mzQC path. Percent transforms below are formatting, not inference.
TMT_METHOD_SOURCES = {
    "spectra_scanned": "levels.<MS>.spectra_scanned",  # centroid only
    "centroid_spectra": "levels.<MS>.centroid_spectra",
    "profile_skipped": "levels.<MS>.profile_skipped",
    "unknown_skipped": "levels.<MS>.unknown_skipped",
    "profile_spectra_scanned": "levels.<MS>.profile_spectra_scanned",
    "profile_spectra_evaluable": "levels.<MS>.profile_spectra_evaluable",
    "profile_spectra_without_reporter_region": "levels.<MS>.profile_spectra_without_reporter_region",
    "profile_spectra_invalid_samples": "levels.<MS>.profile_spectra_invalid_samples",
    "profile_peak_picker_failures": "levels.<MS>.profile_peak_picker_failures",
    "profile_spectra_with_no_picked_peaks": "levels.<MS>.profile_spectra_with_no_picked_peaks",
    "non_profile_spectra_excluded": "levels.<MS>.non_profile_spectra_excluded",
    "reporter_positive_spectra": "levels.<MS>.methods.<method>.reporter_positive_spectra",
    "reporter_positive_pct": "levels.<MS>.methods.<method>.<positive_fraction> × 100",
}
TMT_CHANNEL_SOURCES = {
    "expected_mz": "levels.<MS>.methods.<method>.channels[].expected_mz",
    "spectra_with_hit": "levels.<MS>.methods.<method>.channels[].spectra_with_hit",
    "eligible_pct": "levels.<MS>.methods.<method>.channels[].<eligible_fraction> × 100",
    "local_sampling_count": (
        "levels.<MS>.methods.<method>.channels[].profile_spectra_with_local_mz_sampling"
    ),
    "local_sampling_pct": (
        "levels.<MS>.methods.<method>.channels[].fraction_of_locally_sampled_profile_spectra × 100"
    ),
    "uncorrected_total_matched_intensity": (
        "levels.<MS>.methods.<method>.channels[].uncorrected_total_matched_intensity"
    ),
    "mean_absolute_mass_error_ppm": (
        "levels.<MS>.methods.<method>.channels[].mean_absolute_mass_error_ppm"
    ),
}

_PROFILE_COUNTS = (
    "profile_spectra_scanned",
    "profile_spectra_evaluable",
    "profile_spectra_without_reporter_region",
    "profile_spectra_invalid_samples",
    "profile_peak_picker_failures",
    "profile_spectra_with_no_picked_peaks",
    "non_profile_spectra_excluded",
)
_CENTROID_COUNTS = (
    "spectra_scanned",
    "centroid_spectra",
    "profile_skipped",
    "unknown_skipped",
)


@dataclass(slots=True)
class TMTReportData:
    method_rows: dict[str, dict[str, Any]] = field(default_factory=dict)
    channel_rows: dict[str, dict[str, Any]] = field(default_factory=dict)
    # One separate heatmap/graph group for each mode, MS level and method.
    # Never aggregate method-catalog counts or centroid/profile denominators.
    percent_plots: dict[tuple[str, str, str], dict[str, dict[str, float | None]]] = (
        field(default_factory=dict)
    )
    invalid_metrics: list[str] = field(default_factory=list)


def _count(value: Any, path: str) -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f"{path} must be a nonnegative integer")
    return value


def _number(value: Any, path: str, *, nullable: bool = False) -> float | None:
    # External mzQC field validation uses ValueError for malformed data.
    if value is None and nullable:
        return None
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise ValueError(f"{path} must be a finite number")  # noqa: TRY004
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{path} must be finite")
    return result


def _fraction(value: Any, path: str, *, denominator: int) -> float | None:
    if value is None:
        if denominator > 0:
            raise ValueError(f"{path} missing despite positive denominator")
        return None
    result = _number(value, path)
    if result is None or not 0 <= result <= 1 or denominator == 0:
        raise ValueError(f"{path} inconsistent with its denominator")
    return result * 100


def _read_metric(
    sample_name: str,
    accession: str,
    value: Any,
    result: TMTReportData,
) -> None:
    # Keep structural errors as ValueError: tmt_report_data catches these to
    # record invalid mzQC evidence and leave unaffected runs renderable.
    mode, expected_status = _TMT_SOURCES[accession]
    if not isinstance(value, dict):
        raise ValueError("TMT value must be an object")  # noqa: TRY004
    if (
        type(value.get("schema_version")) is not int
        or value["schema_version"] != 1
        or value.get("scientific_status") != expected_status
    ):
        raise ValueError("unsupported TMT schema_version or scientific_status")
    levels = value.get("levels")
    if not isinstance(levels, dict):
        raise ValueError("levels must be an object")  # noqa: TRY004

    # Parse into isolated buffers: invalid input never leaves half-rendered rows.
    method_rows: dict[str, dict[str, Any]] = {}
    channel_rows: dict[str, dict[str, Any]] = {}
    percent_plots: dict[tuple[str, str, str], dict[str, dict[str, float | None]]] = {}
    for ms_level, level in levels.items():
        if ms_level not in ("MS2", "MS3") or not isinstance(level, dict):
            raise ValueError("unknown MS level or malformed level")
        counts = _CENTROID_COUNTS if mode == "centroid" else _PROFILE_COUNTS
        level_counts = {
            key: _count(level.get(key), f"{ms_level}.{key}") for key in counts
        }
        eligible_key = (
            "centroid_spectra" if mode == "centroid" else "profile_spectra_evaluable"
        )
        denominator = level_counts[eligible_key]
        methods = level.get("methods")
        if not isinstance(methods, dict):
            raise ValueError(f"{ms_level}.methods must be an object")  # noqa: TRY004
        for method_name, method in methods.items():
            if (
                not isinstance(method_name, str)
                or not method_name
                or not isinstance(method, dict)
            ):
                raise ValueError("invalid method name or object")
            prefix = f"{sample_name} | {mode} | {ms_level} | {method_name}"
            positive = _count(
                method.get("reporter_positive_spectra"),
                f"{prefix}.reporter_positive_spectra",
            )
            if positive > denominator:
                raise ValueError(f"{prefix}: positives exceed eligible spectra")
            positive_field = (
                "positive_fraction_of_centroid"
                if mode == "centroid"
                else "positive_fraction_of_evaluable_profile"
            )
            positive_pct = _fraction(
                method.get(positive_field),
                f"{prefix}.{positive_field}",
                denominator=denominator,
            )
            if positive_pct is not None and not math.isclose(
                positive_pct / 100, positive / denominator, rel_tol=1e-9, abs_tol=1e-9
            ):
                raise ValueError(
                    f"{prefix}: positive fraction does not match stored count"
                )
            method_rows[prefix] = {
                "source_file": sample_name,
                "mode": mode,
                "ms_level": ms_level,
                "method": method_name,
                **level_counts,
                "reporter_positive_spectra": positive,
                "reporter_positive_pct": positive_pct,
            }
            channels = method.get("channels")
            if not isinstance(channels, list):
                raise ValueError(f"{prefix}.channels must be an array")  # noqa: TRY004
            seen_channels: set[str] = set()
            plot_key = (mode, ms_level, method_name)
            plotted: dict[str, float | None] = {}
            for channel in channels:
                if not isinstance(channel, dict):
                    raise ValueError(f"{prefix}.channels[] must be an object")  # noqa: TRY004
                label = channel.get("channel")
                if not isinstance(label, str) or not label or label in seen_channels:
                    raise ValueError(f"{prefix}: duplicate or invalid channel")
                seen_channels.add(label)
                path = f"{prefix}.{label}"
                expected_mz = _number(channel.get("expected_mz"), f"{path}.expected_mz")
                if expected_mz is None or expected_mz <= 0:
                    raise ValueError(f"{path}: invalid expected m/z")
                hits = _count(
                    channel.get("spectra_with_hit"), f"{path}.spectra_with_hit"
                )
                if hits > denominator:
                    raise ValueError(f"{path}: hits exceed eligible spectra")
                fraction_field = (
                    "fraction_of_centroid_spectra"
                    if mode == "centroid"
                    else "fraction_of_evaluable_profile_spectra"
                )
                hit_pct = _fraction(
                    channel.get(fraction_field),
                    f"{path}.{fraction_field}",
                    denominator=denominator,
                )
                if hit_pct is not None and not math.isclose(
                    hit_pct / 100, hits / denominator, rel_tol=1e-9, abs_tol=1e-9
                ):
                    raise ValueError(
                        f"{path}: hit fraction does not match stored count"
                    )
                local_sampling_count: int | None = None
                local_sampling_pct: float | None = None
                if mode == "profile":
                    local_sampling_count = _count(
                        channel.get("profile_spectra_with_local_mz_sampling"),
                        f"{path}.profile_spectra_with_local_mz_sampling",
                    )
                    if (
                        local_sampling_count > denominator
                        or hits > local_sampling_count
                    ):
                        raise ValueError(f"{path}: invalid local-sampling support")
                    local_sampling_pct = _fraction(
                        channel.get("fraction_of_locally_sampled_profile_spectra"),
                        f"{path}.fraction_of_locally_sampled_profile_spectra",
                        denominator=local_sampling_count,
                    )
                    if local_sampling_pct is not None and not math.isclose(
                        local_sampling_pct / 100,
                        hits / local_sampling_count,
                        rel_tol=1e-9,
                        abs_tol=1e-9,
                    ):
                        raise ValueError(
                            f"{path}: local sampling fraction does not match stored count"
                        )
                intensity = _number(
                    channel.get("uncorrected_total_matched_intensity"),
                    f"{path}.uncorrected_total_matched_intensity",
                )
                if intensity is None or intensity < 0:
                    raise ValueError(f"{path}: invalid intensity")
                ppm = _number(
                    channel.get("mean_absolute_mass_error_ppm"),
                    f"{path}.mean_absolute_mass_error_ppm",
                    nullable=True,
                )
                if ppm is not None and ppm < 0:
                    raise ValueError(f"{path}: invalid mass error")
                if (hits == 0 and ppm is not None) or (hits > 0 and ppm is None):
                    raise ValueError(
                        f"{path}: mass-error support inconsistent with hit count"
                    )
                row_key = f"{prefix} | {label}"
                channel_rows[row_key] = {
                    "source_file": sample_name,
                    "mode": mode,
                    "ms_level": ms_level,
                    "method": method_name,
                    "channel": label,
                    "expected_mz": expected_mz,
                    "spectra_with_hit": hits,
                    "eligible_pct": hit_pct,
                    "local_sampling_count": local_sampling_count,
                    "local_sampling_pct": local_sampling_pct,
                    "uncorrected_total_matched_intensity": intensity,
                    "mean_absolute_mass_error_ppm": ppm,
                }
                plotted[label] = hit_pct
            percent_plots.setdefault(plot_key, {})[sample_name] = plotted
    # Do not overwrite a repeated run/mode/source collision silently.
    if set(method_rows).intersection(result.method_rows) or set(
        channel_rows
    ).intersection(result.channel_rows):
        raise ValueError("duplicate TMT metric for the same run/mode")
    result.method_rows.update(method_rows)
    result.channel_rows.update(channel_rows)
    for key, values in percent_plots.items():
        result.percent_plots.setdefault(key, {}).update(values)


def tmt_report_data(runs: Iterable[Any]) -> TMTReportData:
    """Project only supported prideQC QCPRIDE TMT values into report data.

    Invalid structured records are excluded and identified; all other run/mode
    records remain renderable. No mzQC values are reconstructed from RAW files.
    """
    result = TMTReportData()
    for run in runs:
        for metric in run.metrics:
            if metric.accession not in _TMT_SOURCES:
                continue
            try:
                _read_metric(run.sample_name, metric.accession, metric.value, result)
            except ValueError as exc:
                result.invalid_metrics.append(
                    f"{run.sample_name} {metric.accession}: {exc}"
                )
    return result


def tmt_overview_html(report: TMTReportData) -> str:
    """Interactive report-only view of stored run, channel and denominator values."""
    views: dict[tuple[str, str, str], dict[str, Any]] = {}
    index: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    for row in report.method_rows.values():
        key = (row["mode"], row["ms_level"], row["method"])
        run = {
            "name": row["source_file"],
            "eligible": row["centroid_spectra"]
            if row["mode"] == "centroid"
            else row["profile_spectra_evaluable"],
            "positive": row["reporter_positive_spectra"],
            "positive_pct": row["reporter_positive_pct"],
            "diagnostics": {
                k: row[k] for k in (*_CENTROID_COUNTS, *_PROFILE_COUNTS) if k in row
            },
            "channels": {},
        }
        views.setdefault(
            key, {"mode": key[0], "level": key[1], "method": key[2], "runs": []}
        )["runs"].append(run)
        index[(*key, row["source_file"])] = run

    for row in report.channel_rows.values():
        run = index.get(
            (row["mode"], row["ms_level"], row["method"], row["source_file"])
        )
        if run is not None:
            run["channels"][row["channel"]] = {
                "hits": row["spectra_with_hit"],
                "pct": row["eligible_pct"],
                "mz": row["expected_mz"],
                "intensity": row["uncorrected_total_matched_intensity"],
                "ppm": row["mean_absolute_mass_error_ppm"],
                "local_count": row["local_sampling_count"],
                "local_pct": row["local_sampling_pct"],
            }

    # Escape HTML delimiter characters in data before embedding JSON in a script.
    payload = json.dumps(list(views.values()), ensure_ascii=False).replace(
        "<", "\\u003c"
    )
    payload = payload.replace(">", "\\u003e").replace("&", "\\u0026")
    return f"""
<div id="prideqc-tmt-qc" class="prideqc-tmt-qc">
  <p><strong>Exploratory reporter-m/z screening, not corrected TMT quantification.</strong>
  Reference catalogues overlap: selecting a catalogue is a display choice, not
  identification of a TMT kit, channel occupancy or sample identity.</p>
  <div class="prideqc-tmt-controls">
    <label>Reference catalogue <select data-choice="method" aria-label="TMT reference catalogue"></select></label>
    <label>MS level <select data-choice="level" aria-label="MS level"></select></label>
    <label>Peak representation <select data-choice="mode" aria-label="Centroid or profile"></select></label>
  </div>
  <p data-status aria-live="polite"></p>
  <h4>Run summary</h4>
  <div class="prideqc-tmt-scroll"><table class="prideqc-tmt-table" data-summary>
    <thead><tr><th>RAW run</th><th>Eligible scans</th><th>≥3-channel screen positive</th>
    <th>Screen positive / eligible</th></tr></thead><tbody></tbody>
  </table></div>
  <h4>Reporter-compatible peak detection by channel</h4>
  <p>Each cell is the percentage of <strong>eligible scans</strong> with a mass-matched
  reporter peak. Hover for the exact hit count and denominator. 0% is measured;
  N/A means no eligible denominator. A scan may match several channels.</p>
  <div class="prideqc-tmt-scroll"><table class="prideqc-tmt-table prideqc-tmt-matrix" data-matrix>
    <thead></thead><tbody></tbody></table></div>
  <p class="prideqc-tmt-key">Colour: <span>0%</span><span class="prideqc-tmt-mid">50%</span>
  <span class="prideqc-tmt-high">100%</span> &nbsp; Grey: N/A</p>
  <details data-advanced><summary>Advanced evidence for the selected view (counts, exclusions, intensity, mass error)</summary>
    <p>Centroid = native centroid peaks; profile = experimental peak-picked spectra.
    Percentages with no eligible denominator are N/A. Uncorrected signal is not
    corrected TMT quantification. All reference-catalogue and MS-level data are
    retained separately in the complete MultiQC data exports and original mzQC.</p>
    <h5>Scan eligibility / exclusions</h5>
    <div class="prideqc-tmt-scroll"><table class="prideqc-tmt-table" data-diagnostics>
      <thead></thead><tbody></tbody></table></div>
    <h5>Channel-level audit data</h5>
    <div class="prideqc-tmt-scroll"><table class="prideqc-tmt-table" data-channel-detail>
      <thead><tr><th>RAW run</th><th>Channel</th><th>Expected m/z</th><th>Hits</th>
      <th>Hit / eligible</th><th>Locally sampled (profile)</th><th>Hit / locally sampled</th>
      <th>Total matched intensity (uncorrected)</th><th>Mean |ppm error|</th></tr></thead>
      <tbody></tbody></table></div>
  </details>
</div>
<style>
  .prideqc-tmt-qc .prideqc-tmt-controls {{display:flex; flex-wrap:wrap; gap:12px; margin:12px 0}}
  .prideqc-tmt-qc .prideqc-tmt-controls label {{display:flex; flex-direction:column; gap:4px}}
  .prideqc-tmt-qc select {{color:inherit; background:var(--mqc-background-color,#242a30);
    border:1px solid #88949b; border-radius:4px; padding:5px}}
  .prideqc-tmt-scroll {{max-width:100%; overflow-x:auto}}
  .prideqc-tmt-table {{width:100%; border-collapse:collapse; font-size:12px}}
  .prideqc-tmt-table th,.prideqc-tmt-table td {{padding:7px; border-bottom:1px solid #7b848455}}
  .prideqc-tmt-table thead th {{text-align:center; white-space:nowrap}}
  .prideqc-tmt-table td:not(:first-child) {{text-align:center; font-variant-numeric:tabular-nums}}
  .prideqc-tmt-table td:first-child {{max-width:245px; overflow:hidden;
    text-overflow:ellipsis; white-space:nowrap; text-align:left}}
  .prideqc-tmt-matrix td[data-pct] {{min-width:54px; color:#151515; font-weight:600}}
  .prideqc-tmt-matrix td[data-na] {{background:#747f8440; color:inherit}}
  .prideqc-tmt-key {{display:flex; align-items:center; gap:8px; font-size:12px}}
  .prideqc-tmt-key span {{padding:3px 10px; background:rgb(230,245,245); color:#151515}}
  .prideqc-tmt-key .prideqc-tmt-mid {{background:rgb(120,198,177)}}
  .prideqc-tmt-key .prideqc-tmt-high {{background:rgb(31,135,106); color:white}}
</style>
<script type="application/json" id="prideqc-tmt-views">{payload}</script>
<script>
(function () {{
  const root = document.getElementById('prideqc-tmt-qc');
  const data = document.getElementById('prideqc-tmt-views');
  if (!root || !data) return;
  const views = JSON.parse(data.textContent);
  const selectors = Object.fromEntries(['mode','level','method'].map(k =>
    [k, root.querySelector('[data-choice="' + k + '"]')]));
  for (const field of Object.keys(selectors)) {{
    const values = [...new Set(views.map(view => view[field]))];
    for (const val of values.sort()) {{
      const opt = document.createElement('option');
      opt.value = val; opt.textContent = val; selectors[field].appendChild(opt);
    }}
    selectors[field].addEventListener('change', update);
  }}
  // Fixed, labelled presentation preference; never infer kit identity from signals.
  const preferred = views.find(v => v.mode === 'centroid' && v.level === 'MS3'
    && v.method === 'TMT11' && v.runs.some(r => r.eligible > 0)) ||
    views.find(v => v.runs.some(r => r.eligible > 0)) || views[0];
  if (!preferred) return;
  for (const field of Object.keys(selectors)) selectors[field].value = preferred[field];
  const cell = (tag, value, title) => {{
    const el = document.createElement(tag); el.textContent = String(value);
    if (title) el.title = title;
    return el;
  }};
  function update() {{
    const view = views.find(v => Object.keys(selectors).every(k => v[k] === selectors[k].value));
    const sum = root.querySelector('[data-summary] tbody');
    const matrixHead = root.querySelector('[data-matrix] thead');
    const matrixBody = root.querySelector('[data-matrix] tbody');
    const diagnostics = root.querySelector('[data-diagnostics]');
    const auditBody = root.querySelector('[data-channel-detail] tbody');
    sum.replaceChildren(); matrixHead.replaceChildren(); matrixBody.replaceChildren();
    diagnostics.tHead.replaceChildren(); diagnostics.tBodies[0].replaceChildren();
    auditBody.replaceChildren();
    root.querySelector('[data-status]').textContent = view
      ? `${{view.method}} / ${{view.level}} / ${{view.mode}} — display reference only; not a kit identification`
      : 'No mzQC evidence for this combination; choose another view.';
    if (!view) return;
    const channels = [...new Set(view.runs.flatMap(r => Object.keys(r.channels)))];
    const diagnosticKeys = [...new Set(view.runs.flatMap(r => Object.keys(r.diagnostics)))];
    const diagHeader = document.createElement('tr');
    ['RAW run', ...diagnosticKeys].forEach(k => diagHeader.appendChild(cell('th', k)));
    diagnostics.tHead.appendChild(diagHeader);
    const header = document.createElement('tr');
    header.appendChild(cell('th','RAW run'));
    channels.forEach(c => header.appendChild(cell('th', c)));
    matrixHead.appendChild(header);
    view.runs.forEach(run => {{
      const row = document.createElement('tr');
      [run.name, run.eligible, run.positive,
        run.positive_pct === null ? 'N/A' : run.positive_pct.toFixed(2) + '%']
        .forEach((val, i) => row.appendChild(cell('td', val, i === 0 ? run.name : null)));
      sum.appendChild(row);
      const diagRow = document.createElement('tr');
      [run.name, ...diagnosticKeys.map(k => run.diagnostics[k] ?? 'N/A')]
        .forEach(v => diagRow.appendChild(cell('td', v)));
      diagnostics.tBodies[0].appendChild(diagRow);
      const heat = document.createElement('tr');
      heat.appendChild(cell('td', run.name, run.name));
      channels.forEach(channel => {{
        const entry = run.channels[channel];
        const pct = entry ? entry.pct : null;
        const label = pct === null ? 'N/A' : pct.toFixed(1) + '%';
        const extra = entry && entry.local_count !== null
          ? `; locally sampled profiles: ${{entry.local_count}} (${{entry.local_pct === null ? 'N/A' : entry.local_pct.toFixed(1) + '%'}})` : '';
        const tip = `${{run.name}} | ${{channel}}: ${{entry ? entry.hits : 'N/A'}} / ${{run.eligible}} eligible scans (${{label}})${{extra}}`;
        const el = cell('td', label, tip);
        el.setAttribute('aria-label', tip);
        if (pct === null) el.setAttribute('data-na','');
        else {{
          el.setAttribute('data-pct', String(pct));
          const t = pct / 100;
          el.style.backgroundColor = `rgb(${{Math.round(230-199*t)}},${{Math.round(245-110*t)}},${{Math.round(245-139*t)}})`;
        }}
        heat.appendChild(el);
        if (entry) {{
          const detail = document.createElement('tr');
          [run.name, channel, entry.mz, entry.hits, label,
           entry.local_count === null ? 'N/A' : entry.local_count,
           entry.local_pct === null ? 'N/A' : entry.local_pct.toFixed(2) + '%',
           entry.intensity, entry.ppm === null ? 'N/A' : entry.ppm]
            .forEach(v => detail.appendChild(cell('td', v)));
          auditBody.appendChild(detail);
        }}
      }});
      matrixBody.appendChild(heat);
    }});
  }}
  update();
}})();
</script>
"""
