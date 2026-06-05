# FTIR Peak Analyzer

An interactive desktop tool for analysing **FTIR (Fourier-Transform
Infrared)** spectra: paste a spectrum, detect peaks, deconvolve overlapping
bands with a multi-Lorentzian model, and read off the band assignments,
heights, widths, and integrated areas.

> **Author**  Teeranan Nongnual
> **Affiliation**  Department of Chemistry, Faculty of Science,
> Burapha University, Thailand
> **Funding**  Burapha University (BUU) and
> Thailand Science Research and Innovation (TSRI)

---

## Features

- **Paste-in workflow** — copy a two-column (wavenumber, intensity) table
  from Excel/Origin/instrument software and paste with `Ctrl+V`.
- **%T ↔ Absorbance** conversion (`A = 2 − log₁₀(%T)`).
- **Baseline correction** (subtract the data minimum to zero the baseline).
- **Automatic peak detection** with adjustable *min height* and
  *prominence* (SciPy `find_peaks`).
- **Manual peak control** — left-click to add a peak, right-click to
  remove one, optional "snap to nearest local maximum".
- **Multi-Lorentzian deconvolution** of overlapping bands via non-linear
  least squares (`scipy.optimize.curve_fit`), reporting `R²` of the fit.
- **Built-in IR group-assignment table** (O–H, C=O, C–H, fingerprint,
  etc.) so each fitted peak is labelled chemically.
- **Publication-style plot** with smart, non-overlapping peak labels;
  copy the figure to the clipboard or save as 300-dpi PNG.
- **Results table** with position, assignment, height, width, and area —
  one click to copy to the clipboard for further analysis.

## Screens & windows

The UI is a single-window Tkinter app with a coloured ribbon:

1. **DATA & AXIS** (blue) — paste, transform, choose plotting range.
2. **VISUAL CONTROL** (teal) — plot/label font sizes, label spacing.
3. **PEAK DETECTION** (orange) — height/prominence sliders, auto-detect.
4. **MODEL & FIT** (purple) — generate initial Lorentzian guesses,
   run the optimisation, see `R²`.
5. **VISUALIZATION GRAPH** — interactive matplotlib canvas.
6. **ANALYSIS RESULTS TABLE** — fitted peaks with chemical assignments.

## Typical workflow

```
PASTE  →  (%T → ABS)  →  BASE=0  →  AUTO DETECT
                                  ↘  click to add / right-click to remove
       →  GENERATE MODEL  →  START FITTING  →  COPY RESULTS
```

## Requirements

- Windows (the clipboard image copy uses the Win32 clipboard API; the
  rest of the code is portable).
- Python 3.10+ recommended.
- Packages:

  ```
  numpy
  scipy
  matplotlib
  pillow
  ```

  `tkinter` ships with the standard Python installer on Windows.

Install with:

```powershell
python -m venv venv
venv\Scripts\activate
pip install numpy scipy matplotlib pillow
```

## Run

```powershell
python trn-irpeakanalyzer.py
```

A prebuilt Windows executable is also produced from this source
(`ir_peakanalyzer_X.exe`, not tracked in git).

## Data format

Two whitespace- or tab-separated columns, one sample per line:

```
4000.00   98.12
3998.00   97.85
3996.00   97.40
...
```

Either `%T` (transmittance) or absorbance values are accepted; use the
**%T → ABS** button if your data is transmittance.

## Method notes

- **Peak model.** Each band is modelled as a Lorentzian
  `f(x) = A · w² / ((x − x₀)² + w²)`, where `A` is height, `x₀` is centre
  (cm⁻¹), and `w` is the half-width parameter. Integrated area is
  `π · A · w`.
- **Fitting.** All detected peaks are fit simultaneously as a sum of
  Lorentzians using Levenberg–Marquardt (`curve_fit`, `maxfev = 15000`).
- **Goodness of fit.** Reported as `R² = 1 − SSres / SStot`.

## License

Released under a **non-commercial use license** — free to use, study,
modify, and redistribute for academic, educational, or any other
non-commercial purpose. Commercial use requires written permission from
the author. See [`LICENSE`](LICENSE) for the full terms.

## How to cite

If you use FTIR Peak Analyzer in a publication, please cite it as:

> Nongnual, T. (2026). *FTIR Peak Analyzer* (Version 2.01) [Computer
> software]. Department of Chemistry, Faculty of Science, Burapha
> University, Thailand.
> https://github.com/teeranann/FTIR-PeakAnalyzer

BibTeX:

```bibtex
@software{nongnual_ftir_peakanalyzer_2026,
  author       = {Nongnual, Teeranan},
  title        = {FTIR Peak Analyzer},
  version      = {2.01},
  year         = {2026},
  institution  = {Department of Chemistry, Faculty of Science,
                  Burapha University, Thailand},
  note         = {Supported by a research grant from Burapha University
                  and Thailand Science Research and Innovation (TSRI).},
  url          = {https://github.com/teeranann/FTIR-PeakAnalyzer}
}
```

## Acknowledgements

This software was developed with support from a research grant by
**Burapha University (BUU)** and **Thailand Science Research and
Innovation (TSRI)**.
