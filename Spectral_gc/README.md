# BK1 Granger Connectivity Analysis (MATLAB)

This repository contains MATLAB code for spectral Granger causality analysis in the BK1 EEG dataset. The README below documents the Granger connectivity scripts currently being prepared for GitHub from `connectivityProjectCodes/`.

If you only want to share the code shown in the screenshot, the main files are:

- `connectivityProjectCodes/runSaveGrangerData.m`
- `connectivityProjectCodes/saveGrangerData.m`
- `connectivityProjectCodes/extractGrangerBandDataAllSubjectsBK1.m`
- `connectivityProjectCodes/compareGrangerMeditatorsControlsBK1.m`
- `connectivityProjectCodes/compareGrangerTopologyMeditatorsControlsBK1.m`
- `connectivityProjectCodes/runCompareGrangerTopologyBK1.m`
- `connectivityProjectCodes/visualizeGrangerSingleSubjectBK1.m`

## What the code does

These scripts implement a workflow for:

1. loading or generating FieldTrip-ready EEG data for BK1 subjects,
2. computing multivariate spectral Granger causality for each subject and protocol,
3. averaging Granger spectra into frequency bands,
4. summarizing connectivity at both channel and ROI level,
5. comparing meditators and controls with permutation statistics, and
6. visualizing subject-level and group-level connectivity patterns.

The code is written for the local BK1 project structure and is **not fully standalone**. Several scripts assume access to BK1 metadata, montage files, FieldTrip, EEGLAB, and supporting analysis functions in sibling folders.

## Main scripts

| File | Purpose |
| --- | --- |
| `runSaveGrangerData.m` | Batch driver that prepares FieldTrip data if needed and then computes Granger results for all good BK1 subjects and all 8 protocols. |
| `saveGrangerData.m` | Core worker that loads one subject/protocol dataset, splits it into pre/post windows, runs MVAR + Granger analysis in FieldTrip, and saves one output `.mat` file per protocol. |
| `extractGrangerBandDataAllSubjectsBK1.m` | Loads saved Granger spectra across all subjects, averages them within frequency bands, computes channel-level and ROI-level matrices, and saves summary outputs. |
| `compareGrangerMeditatorsControlsBK1.m` | Performs ROI-level group comparison between meditators and controls for pre, post, and post-minus-pre connectivity using permutation tests and FDR correction. |
| `compareGrangerTopologyMeditatorsControlsBK1.m` | Performs graph/topology comparison at the 64-electrode level using node and global metrics such as out-degree, in-degree, asymmetry, hubness, Gini, and hub counts. |
| `runCompareGrangerTopologyBK1.m` | Convenience runner that loops over all 8 protocols and runs both unpaired and paired topology analyses. |
| `visualizeGrangerSingleSubjectBK1.m` | Creates a single-subject figure with pre/post/difference heatmaps, ROI summary, directed pair spectra, and top net-flow electrodes. |

## Default analysis settings

The current code defaults to the following settings:

- Protocols: `EO1`, `EC1`, `G1`, `M1`, `G2`, `EO2`, `EC2`, `M2`
- Bad-eye condition tag: `ep`
- Bad-trial version: `v8`
- Post window: `0.25` to `1.25` seconds
- Matched pre window: approximately `-1.0` to `0` seconds
- MVAR model order: `10`
- MVAR toolbox: `biosig` via FieldTrip
- Frequency range: `0` to `100` Hz in `1` Hz steps
- Default bands:
  - alpha: `7-10` Hz
  - beta: `20-32` Hz
  - highgamma: `30-80` Hz

## Expected folder structure

These scripts infer the project root from their file location, so they expect the broader BK1 workspace layout to stay intact. A typical local setup looks like:

```text
BK1/
├── CommonPrograms/
├── Montages/
├── ProjectDhyaanBK1Programs/
│   ├── commonAnalysisCodes/
│   ├── connectivityProjectCodes/
│   └── ...
├── eeglab/
├── fieldtrip-20260211/
└── data/
    └── ftData/
```

Important dependencies used by the scripts:

- MATLAB
- FieldTrip (`fieldtrip-20260211`)
- EEGLAB topoplot support for topology figures
- BK1 metadata in `commonAnalysisCodes/informationFiles/`
- Montage labels in `Montages/Layouts/actiCap64_UOL/`

## Outputs

The analysis writes derived files under `connectivityProjectCodes/`, mainly:

- `savedDataGranger/<subject>/<protocol>_ep_v8_granger.mat`
- `savedDataGranger/summary/`
- `savedDataGranger/groupStats/`
- `savedDataGranger/topologyStats/`
- `savedDataGranger/figures/singleSubject/`

These are generated results rather than source code, so they are typically not needed in a code-only GitHub upload.

## Recommended workflow

### 1. Generate or refresh per-subject Granger outputs

```matlab
cd connectivityProjectCodes
runSaveGrangerData(0)
```

Set `runSaveGrangerData(1)` if you also want to regenerate the intermediate FieldTrip data first.

### 2. Create a summary for one protocol

```matlab
results = extractGrangerBandDataAllSubjectsBK1('M1');
```

### 3. Compare meditators vs controls

ROI-level comparison:

```matlab
statsROI = compareGrangerMeditatorsControlsBK1('M1');
```

Topology-level comparison:

```matlab
statsTopo = compareGrangerTopologyMeditatorsControlsBK1('M1');
```

Run all protocols:

```matlab
runCompareGrangerTopologyBK1
```

### 4. Visualize a single subject

```matlab
visualizeGrangerSingleSubjectBK1('028HB','M1','alpha');
```

## Notes for GitHub upload

- The code in this folder depends on the rest of the BK1 workspace, so uploading only these `.m` files is best for code sharing, but not enough for full reproducibility by itself.
- If you want a clean code repository, upload the scripts and this `README.md`, but avoid committing large derived `.mat` outputs and figures.
- A minimal `.gitignore` is included to keep generated result folders out of Git by default.
