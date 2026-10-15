"""mzQC -> TMT reporter rendering contract, without RAW or sidecar inputs."""

from __future__ import annotations

import copy
import json
import re
from pathlib import Path
from types import SimpleNamespace

from pmultiqc.modules.mzqc import parse_mzqc_document
from pmultiqc.modules.mzqc.tmt_reporting import (
    CENTROID_ACCESSION,
    PROFILE_ACCESSION,
    TMT_CHANNEL_SOURCES,
    TMT_METHOD_SOURCES,
    tmt_overview_html,
    tmt_report_data,
)

FIXTURE = Path(__file__).parent / "resources" / "mzqc" / "prideqc_tmt_v1.mzQC"


def _check(condition: bool) -> None:
    if not condition:
        raise AssertionError("TMT mzQC reporting parity violation")


def _fixture():
    runs = parse_mzqc_document(FIXTURE)
    _check(len(runs) == 1)
    return runs[0]


def _metric(run, accession):
    return next(m for m in run.metrics if m.accession == accession)


def test_tmt_run_overview_keeps_structured_evidence_out_of_table_cells():
    """An actual TMT report must not print full nested mzQC JSON in Run Overview."""
    from pmultiqc.modules.mzqc import mzqc as module

    run = _fixture()
    overview = module.MzQCModule._run_overview_data(SimpleNamespace(runs=[run]))
    label = overview[run.sample_name]["provenance"]
    _check(label == "TMT reporter evidence: see TMT QC sections")
    _check("schema_version" not in label)
    _check("levels" not in label)
    _check(len(label) < 100)
    # The original mzQC evidence and dedicated TMT projections are not altered.
    _check(_metric(run, CENTROID_ACCESSION).value["levels"])
    _check(_metric(run, PROFILE_ACCESSION).value["levels"])
    _check(bool(tmt_report_data([run]).channel_rows))


def test_overview_retains_scalar_provenance_and_summarizes_other_structures():
    from pmultiqc.modules.mzqc.mzqc import MzQCMetric, MzQCRun

    run = MzQCRun(
        sample_name="test.raw",
        source_path=Path("test.raw.mzQC"),
        metadata={},
        metrics=(
            MzQCMetric("QCPRIDE:scalar", "acquisition method provenance", "inferred"),
            MzQCMetric(
                "QCPRIDE:nested", "other evidence", {"very_long_key": [1, 2, 3]}
            ),
        ),
    )
    label = run.provenance
    _check(label is not None)
    _check('acquisition method provenance: "inferred"' in label)
    _check("other evidence: structured value (see source mzQC)" in label)
    _check("very_long_key" not in label)


def test_tmt_all_exported_numbers_equal_mzqc_source_fields():
    run = _fixture()
    report = tmt_report_data([run])
    _check(not report.invalid_metrics)
    for accession, mode in (
        (CENTROID_ACCESSION, "centroid"),
        (PROFILE_ACCESSION, "profile"),
    ):
        source = _metric(run, accession).value
        for ms_level, level in source["levels"].items():
            for method, value in level["methods"].items():
                prefix = f"{run.sample_name} | {mode} | {ms_level} | {method}"
                row = report.method_rows[prefix]
                for key in (
                    "spectra_scanned",
                    "centroid_spectra",
                    "profile_skipped",
                    "unknown_skipped",
                    "profile_spectra_scanned",
                    "profile_spectra_evaluable",
                    "profile_spectra_without_reporter_region",
                    "profile_spectra_invalid_samples",
                    "profile_peak_picker_failures",
                    "profile_spectra_with_no_picked_peaks",
                    "non_profile_spectra_excluded",
                ):
                    if key in level:
                        _check(row[key] == level[key])
                _check(
                    row["reporter_positive_spectra"]
                    == value["reporter_positive_spectra"]
                )
                method_fraction_key = (
                    "positive_fraction_of_centroid"
                    if mode == "centroid"
                    else "positive_fraction_of_evaluable_profile"
                )
                fraction = value[method_fraction_key]
                _check(
                    row["reporter_positive_pct"]
                    == (None if fraction is None else 100 * fraction)
                )
                for channel in value["channels"]:
                    channel_row = report.channel_rows[
                        f"{prefix} | {channel['channel']}"
                    ]
                    for key in (
                        "expected_mz",
                        "spectra_with_hit",
                        "uncorrected_total_matched_intensity",
                        "mean_absolute_mass_error_ppm",
                    ):
                        _check(channel_row[key] == channel[key])
                    eligible_field = (
                        "fraction_of_centroid_spectra"
                        if mode == "centroid"
                        else "fraction_of_evaluable_profile_spectra"
                    )
                    eligible = channel[eligible_field]
                    _check(
                        channel_row["eligible_pct"]
                        == (None if eligible is None else 100 * eligible)
                    )
                    if mode == "profile":
                        _check(
                            channel_row["local_sampling_count"]
                            == channel["profile_spectra_with_local_mz_sampling"]
                        )
                        local_fraction = channel[
                            "fraction_of_locally_sampled_profile_spectra"
                        ]
                        _check(
                            channel_row["local_sampling_pct"]
                            == (
                                None if local_fraction is None else 100 * local_fraction
                            )
                        )
                    else:
                        _check(channel_row["local_sampling_count"] is None)
                        _check(channel_row["local_sampling_pct"] is None)


def test_tmt_fraction_denominators_missing_remain_na_but_observed_zero_is_zero():
    report = tmt_report_data([_fixture()])
    centroid_ms3 = next(
        v for k, v in report.channel_rows.items() if "centroid | MS3" in k
    )
    zero_hits = next(
        v for k, v in report.channel_rows.items() if "centroid | MS2 | TMT6 | 127" in k
    )
    missing_profile = next(
        v for k, v in report.channel_rows.items() if "profile | MS2 | TMT6 | 127" in k
    )
    _check(centroid_ms3["eligible_pct"] is None)
    _check(zero_hits["eligible_pct"] == 0)
    _check(missing_profile["local_sampling_pct"] is None)
    _check(missing_profile["eligible_pct"] == 0)
    # A missing denominator never becomes a charted 0% bar.
    _check(
        report.percent_plots[("centroid", "MS3", "TMT6")]["prideqc-tmt-demo.raw"]["126"]
        is None
    )


def test_tmt_catalogs_ms_levels_modes_remain_independent():
    report = tmt_report_data([_fixture()])
    _check(("centroid", "MS2", "TMT6") in report.percent_plots)
    _check(("centroid", "MS2", "TMT11") in report.percent_plots)
    _check(("profile", "MS2", "TMT6") in report.percent_plots)
    _check(len(report.channel_rows) == 7)
    _check(
        all("total_unique_spectra" not in row for row in report.method_rows.values())
    )


def test_tmt_rejects_bad_or_unexpected_metric_without_partial_output():
    original = _fixture()
    metrics = list(original.metrics)
    centroid = copy.deepcopy(_metric(original, CENTROID_ACCESSION))
    # MzQCMetric is frozen; make an editable metric adapter for this test.
    evidence = copy.deepcopy(centroid.value)
    evidence["levels"]["MS2"]["methods"]["TMT6"]["channels"][0]["spectra_with_hit"] = 50
    bad = SimpleNamespace(accession=CENTROID_ACCESSION, value=evidence)
    valid = _metric(original, PROFILE_ACCESSION)
    report = tmt_report_data([SimpleNamespace(sample_name="bad", metrics=[bad, valid])])
    _check(len(report.invalid_metrics) == 1)
    _check(all("centroid" not in name for name in report.method_rows))
    _check(any("profile" in name for name in report.method_rows))
    _check(len(metrics) == 2)


def test_tmt_only_supported_accessions_and_schema_are_rendered():
    original = _fixture()
    unknown = SimpleNamespace(
        accession="QCPRIDE:OTHER", value=_metric(original, CENTROID_ACCESSION).value
    )
    wrong_version = copy.deepcopy(_metric(original, CENTROID_ACCESSION).value)
    wrong_version["schema_version"] = 2
    invalid = SimpleNamespace(accession=CENTROID_ACCESSION, value=wrong_version)
    report = tmt_report_data(
        [SimpleNamespace(sample_name="unknown", metrics=[unknown, invalid])]
    )
    _check(not report.method_rows)
    _check(len(report.invalid_metrics) == 1)


def test_tmt_invalid_field_types_remain_fail_closed_value_errors():
    """Bad mzQC JSON value types must not escape the reporting loader."""
    original = _fixture()
    good = _metric(original, PROFILE_ACCESSION)
    for invalid_value in (
        [],
        {
            "schema_version": 1,
            "scientific_status": "EXPLORATORY_NOT_FOR_AUTOMATIC_SDRF_WRITEBACK",
            "levels": [],
        },
    ):
        bad = SimpleNamespace(accession=CENTROID_ACCESSION, value=invalid_value)
        report = tmt_report_data(
            [SimpleNamespace(sample_name="invalid-types", metrics=[bad, good])]
        )
        _check(len(report.invalid_metrics) == 1)
        _check(bool(report.channel_rows))
        _check(all("centroid" not in key for key in report.method_rows))


def test_tmt_report_numeric_columns_have_explicit_mzqc_paths():
    method_rows = tmt_report_data([_fixture()]).method_rows
    channel_rows = tmt_report_data([_fixture()]).channel_rows
    method_numeric = {
        k
        for row in method_rows.values()
        for k, v in row.items()
        if isinstance(v, (int, float)) and not isinstance(v, bool)
    }
    channel_numeric = {
        k
        for row in channel_rows.values()
        for k, v in row.items()
        if isinstance(v, (int, float)) and not isinstance(v, bool)
    }
    _check(method_numeric == set(TMT_METHOD_SOURCES))
    _check(channel_numeric == set(TMT_CHANNEL_SOURCES))
    _check(all("levels." in path for path in TMT_METHOD_SOURCES.values()))
    _check(all("channels[]" in path for path in TMT_CHANNEL_SOURCES.values()))


def test_tmt_inconsistent_count_and_fraction_abstains():
    run = _fixture()
    evidence = copy.deepcopy(_metric(run, CENTROID_ACCESSION).value)
    evidence["levels"]["MS2"]["methods"]["TMT6"]["channels"][0][
        "fraction_of_centroid_spectra"
    ] = 0.1
    bad = SimpleNamespace(accession=CENTROID_ACCESSION, value=evidence)
    report = tmt_report_data([SimpleNamespace(sample_name="mismatch", metrics=[bad])])
    _check(not report.method_rows)
    _check(len(report.invalid_metrics) == 1)


def test_tmt_renderer_shows_source_values_and_omits_na_plots(monkeypatch):
    from pmultiqc.modules.mzqc import mzqc as module

    sections = []
    outputs = {}
    plots = []

    def fake_plot(data, **kwargs):
        plots.append((data, kwargs["pconfig"]["id"]))
        return "rendered"

    monkeypatch.setattr(module.table, "plot", fake_plot)
    monkeypatch.setattr(module.bargraph, "plot", fake_plot)
    fake = SimpleNamespace(
        runs=[_fixture()],
        log=SimpleNamespace(warning=lambda *args: None),
        add_section=lambda **kwargs: sections.append(kwargs),
        write_data_file=lambda data, name: outputs.update({name: data}),
    )
    module.MzQCModule._add_tmt_sections(fake)
    _check(any(s["name"].startswith("TMT QC —") for s in sections))
    _check("multiqc_mzqc_tmt_source_map" in outputs)
    _check("multiqc_mzqc_tmt_channels" in outputs)
    _check(not plots)  # No duplicated tables or catalogue-by-catalogue plots.
    _check(not any("channel_pct" in pid for _, pid in plots))
    overview = next(s["content"] for s in sections if s["name"].startswith("TMT QC"))
    _check('data-choice="method"' in overview)
    _check('data-choice="mode"' in overview)
    _check('data-choice="level"' in overview)
    _check("data-matrix" in overview)
    _check("data-advanced" in overview)
    _check("data-channel-detail" in overview)
    _check("data-diagnostics" in overview)
    _check("No mzQC evidence for this combination" in overview)


def test_tmt_overview_data_is_strictly_source_mapped_and_null_safe():
    report = tmt_report_data([_fixture()])
    html = tmt_overview_html(report)
    match = re.search(
        r'<script type="application/json" id="prideqc-tmt-views">(.*?)</script>',
        html,
        re.DOTALL,
    )
    _check(match is not None)
    views = json.loads(match.group(1))
    centroid_ms2 = next(
        v
        for v in views
        if (v["mode"], v["level"], v["method"]) == ("centroid", "MS2", "TMT6")
    )
    run = centroid_ms2["runs"][0]
    source = report.method_rows["prideqc-tmt-demo.raw | centroid | MS2 | TMT6"]
    _check(run["eligible"] == source["centroid_spectra"])
    _check(run["positive"] == source["reporter_positive_spectra"])
    _check(run["positive_pct"] == source["reporter_positive_pct"])
    _check(run["channels"]["126"]["hits"] == 9)
    _check(run["channels"]["126"]["pct"] == 90.0)
    _check(
        run["channels"]["126"]["mz"]
        == report.channel_rows["prideqc-tmt-demo.raw | centroid | MS2 | TMT6 | 126"][
            "expected_mz"
        ]
    )
    _check(run["channels"]["127"]["pct"] == 0.0)
    _check(run["diagnostics"]["centroid_spectra"] == source["centroid_spectra"])
    unavailable = next(
        v
        for v in views
        if (v["mode"], v["level"], v["method"]) == ("centroid", "MS3", "TMT6")
    )
    _check(unavailable["runs"][0]["channels"]["126"]["pct"] is None)
    # The default view choice is a labelled preference, not a data-derived kit call.
    _check("display reference only; not a kit identification" in html)
    _check("entry ? entry.hits" in html)
    _check("run.eligible" in html)


def test_tmt_overview_escapes_html_breakout_in_input_name():
    report = tmt_report_data([_fixture()])
    original = next(iter(report.method_rows.values()))["source_file"]
    poison = "</script><script>alert(1)</script>"
    for row in report.method_rows.values():
        if row["source_file"] == original:
            row["source_file"] = poison
    html = tmt_overview_html(report)
    _check(poison not in html)
    _check(r"\u003c/script\u003e" in html)


def test_tmt_and_standard_mzqc_metrics_render_together(monkeypatch):
    """Adding TMT records must not suppress the existing scalar QC sections."""
    from pmultiqc.modules.mzqc import mzqc as module

    path = FIXTURE.with_name("prideqc_tmt_with_standard_qc_v1.mzQC")
    runs = parse_mzqc_document(path)
    _check(len(runs) == 1)
    run = runs[0]
    # This is an entirely synthetic fixture; the numbers exercise rendering,
    # not scientific validity of a real acquisition.
    _check(len(run.metrics) > len(_fixture().metrics))
    _check(
        next(m.value for m in run.metrics if m.name == "number of MS2 spectra") == 12
    )
    sections = []
    plots = []
    outputs = {}

    def fake_plot(data, **kwargs):
        plots.append((data, kwargs["pconfig"]["id"]))
        return "rendered"

    monkeypatch.setattr(module.table, "plot", fake_plot)
    monkeypatch.setattr(module.bargraph, "plot", fake_plot)
    fake = SimpleNamespace(
        runs=runs,
        log=SimpleNamespace(warning=lambda *args: None),
        add_section=lambda **kwargs: sections.append(kwargs),
        write_data_file=lambda data, name: outputs.update({name: data}),
    )
    for title, anchor, categories in (
        ("Chromatography", "mzqc-chromatography", {"chromatography"}),
        ("MS1 / MS2 Summary", "mzqc-spectra", {"spectra"}),
        ("Acquisition", "mzqc-acquisition", {"acquisition"}),
        ("Other Mass Accuracy / Precision Metrics", "mzqc-mass-other", {"mass"}),
    ):
        module.MzQCModule._add_category_table(fake, title, anchor, categories)
    module.MzQCModule._add_tmt_sections(fake)
    names = {section["name"] for section in sections}
    _check("Chromatography" in names)
    _check("MS1 / MS2 Summary" in names)
    _check("Acquisition" in names)
    _check("Other Mass Accuracy / Precision Metrics" in names)
    _check(any(name.startswith("TMT QC —") for name in names))
    _check(not any("channel_pct" in plot_id for _, plot_id in plots))
    _check("multiqc_mzqc_tmt_methods" in outputs)
    _check("multiqc_mzqc_tmt_channels" in outputs)
    # Counts and reporter percentages remain tied to the same mzQC run.
    _check(
        outputs["multiqc_mzqc_tmt_channels"][
            "prideqc-tmt-demo.raw | centroid | MS2 | TMT6 | 126"
        ]["spectra_with_hit"]
        == 9
    )
