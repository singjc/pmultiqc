# prideQC TMT reporter-ion screening in mzQC reports

The existing `mzqc` module in pmultiqc can display prideQC's optional TMT reporter-ion screening evidence. No RAW file, peak detection, identification or `.summary.json` sidecar is required at report time:

```bash
multiqc --module mzqc --outdir report results/
```

The module reads two structured, local-vocabulary mzQC quality metrics, when present:

| mzQC accession | Record | Scientific scope |
| --- | --- | --- |
| `QCPRIDE:D9D5612FDC49DA507254` | TMT reporter-ion channel evidence | Native centroid MS2/MS3 |
| `QCPRIDE:94464A800A02DA828D7E` | TMT experimental profile reporter-ion evidence | Ephemerally peak-picked profile MS2/MS3 |

The TMT section is now **one interactive overview**, not a separate plot for every reference catalogue. Choose the reference catalogue (TMT6, TMT10, TMT11, TMTpro16 or TMTpro18), MS level (MS2 or MS3), and spectrum representation (centroid or experimental profile). When an eligible **TMT11 / MS3 / centroid** view exists, it is the initial display *by reporting preference only* — not a discovered kit identity. Otherwise, the initial view is the first combination with eligible spectra. Combinations absent from mzQC are reported as unavailable, never fabricated.

The selected view shows a per-run summary: eligible scans, scans passing the exploratory three-channel screen and the **stored reporter-positive fraction × 100**. The channel-by-run heatmap shows the **stored per-channel hit fraction × 100**; hovering shows each hit count, denominator and (for profile spectra) the local sampling count. **0% is not N/A**. The details disclosure shows the selected view's exclusions, expected m/z, hit counts, local sampling, uncorrected matched intensity and mean absolute ppm mass error. The full, nonpooled data from **all** reference catalogues and modes remain in the original mzQC files and in `multiqc_mzqc_tmt_methods`, `multiqc_mzqc_tmt_channels` and `multiqc_mzqc_tmt_source_map` report data exports. No separate channel charts or permanently expanded thousands-row tables are generated.

All numerical data originate from the mzQC `value.levels` structures, and percentage conversion is display-only (`stored fraction × 100`). No new reporter hits, thresholds or biological interpretations are computed by pmultiqc.

**Denominators matter:** Centroid percentages use `centroid_spectra`, profile percentages use `profile_spectra_evaluable`, and locally sampled profile percentages use each channel's `profile_spectra_with_local_mz_sampling`. If a denominator is zero the fraction is unavailable (`N/A`), not zero. MS2 and MS3, centroid and profile, and multiple overlapping TMT reference catalogues are never pooled.

**Interpretation:** These metrics are exploratory, identification-free peak screening, *not* TMT isotope-impurity-corrected quantification. Catalogue matches overlap and do not identify the kit, active sample channels, sample assignment, or SDRF labels. The reported number of reporter-positive spectra reflects the threshold stored in prideQC's original mzQC evidence; pmultiqc does not reclassify scans.

Malformed or unsupported local TMT metric records are skipped with a warning rather than being converted into apparently valid QC figures; ordinary mzQC reporting remains operational.

## Synthetic smoke fixtures versus a complete prideQC mzQC report

`tests/resources/mzqc/prideqc_tmt_v1.mzQC` contains **only the two
structured TMT evidence metrics**. It is intentionally minimal for edge-case
tests, so it does **not** produce standard Chromatography, MS1 / MS2 Summary,
Acquisition or mass-accuracy tables. Their absence is expected, not evidence
that the TMT extension hides standard QC sections.

To exercise **standard QC and TMT together**, run the separate, completely
synthetic combined fixture:

```bash
multiqc --module mzqc --force --outdir /tmp/pmultiqc-tmt-combined \
  tests/resources/mzqc/prideqc_tmt_with_standard_qc_v1.mzQC
```

The combined fixture adds illustrative mzQC scalar values (chromatography
duration, MS1/MS2 scan counts, acquisition mode, MS2 isolation width, m/z range, precursor
precision) without modifying the TMT evidence. These values are invented
**test data**, not measurements or benchmarks, and must not be represented as
results from a real RAW file.

For actual DDA/DIA/TMT reports, supply **complete prideQC-generated mzQC
files** to `multiqc --module mzqc`. The renderer only displays existing
metrics; it does not infer absent values from the TMT annotation or RAW data.
Other independent pmultiqc modules may require their own input formats and
are not automatically activated by `--module mzqc`.
