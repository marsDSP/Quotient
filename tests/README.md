# Quotient filter tests

Numeric verification of `MarsDSP::Filters::TwoPoleSVF`
(`source/dsp/filters/StateVariable.h`) against Andrew Simper / Cytomic,
*"Solving the continuous SVF equations using trapezoidal integration and
equivalent currents"* (`docs/papers/SvfLinearTrapOptimised2.pdf`).

## Layout

- `numeric/` — Catch2 v3 numeric tests (`SvfNumericTests.cpp`) plus an
  independent, paper-derived reference (`SvfReference.h`) used as the oracle.
  Also implements the chart-data exporter (`SvfChartData.h`).
- `harnesses/` — the filter test harness (`FilterTestHarness.cpp`). It owns
  Catch2's `main`, runs the numeric suite, and on `--export-charts` writes CSVs.
- `charts/` — `plot_svf.py` renders PNG charts from the exported CSVs (written
  to `charts/data/`).

## Build & run (standalone — recommended, no JUCE/KFR)

```sh
cmake -S tests -B build-filter-test
cmake --build build-filter-test
ctest --test-dir build-filter-test --output-on-failure
```

Generate the charts (runs the suite, exports CSVs, renders PNGs into `charts/`):

```sh
cmake --build build-filter-test --target charts
```

Or manually:

```sh
./build-filter-test/harnesses/quotient_filter_tests --export-charts
python3 tests/charts/plot_svf.py
```

## Build as part of the main project

```sh
cmake -S . -B build -DQUOTIENT_BUILD_TESTS=ON
cmake --build build --target quotient_filter_tests
```

## What is verified

1. `magnitude()` matches the paper's exact discrete transfer function.
2. `processSample()` matches an independent transcription of the paper's
   bounded trapezoidal algorithm, sample-for-sample.
3. A black-box impulse-response DFT matches the analytic response.
4. Characteristic-frequency behaviour, stability, `OnePoleTPT` and `TiltShelf`
   sanity (the latter two are not part of the paper).

### Known, intentional deviation

`BandPass` is normalised to a **0 dB peak** (`m1 = 1/Q`). The paper's generic
`band` response uses `m1 = 1`, i.e. a peak gain of `Q`. The two differ by
exactly a factor of `Q`; this is asserted explicitly in
`"BandPass is normalised to 0 dB peak; paper band peaks at Q"`.
