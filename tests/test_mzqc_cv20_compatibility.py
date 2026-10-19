"""Regression checks for the prideQC 20-metric mzQC contract."""

from pathlib import Path

from pmultiqc.modules.mzqc.mzqc import (
    MzQCMetric,
    MzQCModule,
    MzQCRun,
    _experiment_group_features,
    _metric_data,
    _run_duration_seconds,
    _run_group_assignments,
    classify_metric,
)


def _check(condition: bool) -> None:
    if not condition:
        raise AssertionError("CV20 mzQC compatibility regression")


def _run(name: str, *metrics: MzQCMetric, metadata: dict | None = None) -> MzQCRun:
    return MzQCRun(name, Path(f"{name}.mzQC"), metadata or {}, tuple(metrics))


def _duration(value: float, accession: str = "MS:4000067") -> MzQCMetric:
    name = "MS run duration" if accession == "MS:4000067" else "chromatography duration"
    return MzQCMetric(accession, name, value, "UO:0000010", "second")


def test_cv20_duration_read_from_official_accession_and_historical_name() -> None:
    current = _run("current", _duration(5400))
    historical = _run("historical", _duration(7200, "MS:4000053"))
    unrelated = _run("unrelated", MzQCMetric("MS:4000053", "unrelated measurement", 60))

    _check(_run_duration_seconds(current) == 5400)
    _check(_run_duration_seconds(historical) == 7200)
    _check(_run_duration_seconds(unrelated) is None)

    _, features = _experiment_group_features([current, historical, unrelated])
    _check(features["current"]["rt"] == 5400)
    _check(features["historical"]["rt"] == 7200)
    _check(features["unrelated"]["rt"] is None)


def test_retired_isolation_is_not_used_for_experiment_group_inference() -> None:
    runs = [
        _run(f"a{i}", _duration(5400), MzQCMetric("QCPRIDE:old", "IsolationWidth_MS2_Median", 1))
        for i in range(4)
    ] + [
        _run(f"b{i}", _duration(5400), MzQCMetric("QCPRIDE:old", "IsolationWidth_MS2_Median", 50))
        for i in range(4)
    ]
    _, features = _experiment_group_features(runs)
    _check(all("isolation" not in values for values in features.values()))
    groups, summaries, evidence = _run_group_assignments(runs)
    _check(len(summaries) == 1)
    _check(set(groups.values()) == {"Experiment group 1"})
    _check(not any("isolation width" in item for item in evidence))
    _check("isolation_median_th" not in summaries["Experiment group 1"])


def test_official_dia_and_observed_bounds_display_exact_mzqc_values() -> None:
    run = _run(
        "dia",
        MzQCMetric("MS:4000194", "DIA isolation window count", 17),
        MzQCMetric("MS:4000195", "DIA isolation window m/z widths", [8, 26]),
        MzQCMetric("MS:4000222", "MS1 observed m/z bounds", [352, 1698]),
    )
    data, headers = _metric_data([run])
    _check(data["dia"]["MS_4000194"] == 17)
    _check(data["dia"]["MS_4000195_min"] == 8)
    _check(data["dia"]["MS_4000195_max"] == 26)
    _check(data["dia"]["MS_4000222_min"] == 352)
    _check(data["dia"]["MS_4000222_max"] == 1698)
    _check(headers["MS_4000195_min"]["title"] == "DIA isolation window width — minimum")
    _check(headers["MS_4000195_max"]["title"] == "DIA isolation window width — maximum")


def test_faims_observed_metadata_does_not_generate_new_qc() -> None:
    metadata = {
        "inputFiles": [
            {"name": "run.raw", "fileProperties": [
                {"accession": "MS:1001581", "name": "FAIMS compensation voltage", "value": -45.0,
                 "unit": {"accession": "UO:0000218", "name": "volt"}},
                {"accession": "MS:1001581", "name": "FAIMS compensation voltage", "value": -65.0,
                 "unit": {"accession": "UO:0000218", "name": "volt"}},
                {"accession": "MS:1003450", "name": "FAIMS ramp start", "value": -35},
                {"accession": "MS:1001581", "name": "bad unit", "value": -2,
                 "unit": {"accession": "UO:0000010", "name": "second"}},
            ]}
        ]
    }
    run = _run("faims", _duration(3600), metadata=metadata)
    module = object.__new__(MzQCModule)
    module.runs = [run]
    overview = module._run_overview_data()["faims"]
    _check(overview["faims_cv"] == "-45.0 V; -65.0 V")
    qc, _ = _metric_data([run])
    _check(qc["faims"] == {"MS_4000067": 3600})


def test_legacy_duration_is_not_labeled_as_new_metric() -> None:
    run = _run("old", _duration(4000, "MS:4000053"))
    data, headers = _metric_data([run])
    _check(data["old"] == {"MS_4000053": 4000})
    _check("MS_4000067" not in headers)


def test_cv20_metric_categories_follow_official_semantics() -> None:
    _check(classify_metric(_duration(4800)) == "chromatography")
    _check(
        classify_metric(MzQCMetric("MS:4000192", "MS1 median cycle time", 1.2))
        == "acquisition"
    )
