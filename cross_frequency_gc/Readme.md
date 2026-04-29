# Cross-Frequency Granger Causality

## Overview

This folder contains a MATLAB workflow for BK1 cross-frequency phase-to-amplitude Granger causality (CF-GC). The code computes subject-level CF-GC matrices, saves them to disk, and then provides visualization, matched-pair comparison, protocol-difference, secondary summary, and net-flow analysis layers on top of those saved outputs.

The workflow is designed around BK1 data conventions and expects the surrounding project folders to be available, especially:

- `ProjectDhyaanBK1Programs`
- `CommonPrograms`
- `Montages`
- `fieldtrip-20260211`

## Expected Inputs And Outputs

Inputs:

- FieldTrip-ready subject files under the folder used by `computeCrossFrequencyPhaseAmplitudeGrangerBK1`.
- BK1 subject metadata and matched-pair definitions from `ProjectDhyaanBK1Programs/commonAnalysisCodes/informationFiles`.

Outputs:

- `savedDataCrossFreqGranger/` for subject-level CF-GC MAT files
- `analysis/` for visualization and comparison outputs
- `analysis_2/` for secondary summary MAT files
- `analysis_net_flow/` for net-flow outputs

## Typical Usage

1. Generate subject-level CF-GC files:

```matlab
runSaveCrossFrequencyPhaseAmplitudeGrangerBK1('all',{'pre','post'},'ep','v8');
```

2. Run the main visualization/comparison launcher:

```matlab
runAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1;
```

3. Or call a specific entry point directly:

```matlab
visualizeCrossFrequencyPhaseAmplitudeGrangerSingleSubjectBK1('019CKa','M1','ep','v8');
compareCrossFrequencyPhaseAmplitudeGrangerMultiplePairsBK1('all',{'M1','M2'},'ep','v8');
compareCrossFrequencyPhaseAmplitudeGrangerProtocolDifferenceBK1('G2','G1','all','ep','v8');
```

4. Optional follow-up analyses:

```matlab
runAnalysis2CrossFrequencyPhaseAmplitudeGrangerBK1('all','all','ep','v8');
runNetFlowAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1('all','all','ep','v8');
```

## Notes

- The core computation uses a linearized phase representation by default, so the CF-GC results should be treated as exploratory unless you add stronger confirmatory validation.
- The comparison workflow is built around matched meditator-control pairs.
- Most analysis functions support cached outputs and will reuse saved summaries unless `forceRebuild` is enabled.

## File Guide

### Main Entry Points

- `computeCrossFrequencyPhaseAmplitudeGrangerBK1.m`: Computes one subject/protocol/condition CF-GC result from FieldTrip input.
- `runSaveCrossFrequencyPhaseAmplitudeGrangerBK1.m`: Batch-runs the subject-level computation and saves one MAT file per subject/protocol/condition.
- `runAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1.m`: Interactive launcher that lets you switch between visualization and comparison workflows by editing one configuration block.
- `visualizeCrossFrequencyPhaseAmplitudeGrangerSingleSubjectBK1.m`: Creates descriptive PRE, POST, and POST-PRE figures for one subject.
- `visualizeCrossFrequencyPhaseAmplitudeGrangerMultipleSubjectsBK1.m`: Repeats the single-subject visualization workflow for multiple validated subjects and saves a combined summary.
- `compareCrossFrequencyPhaseAmplitudeGrangerSinglePairBK1.m`: Builds descriptive figures for one matched meditator-control pair without group-level inference.
- `compareCrossFrequencyPhaseAmplitudeGrangerMultiplePairsBK1.m`: Runs the main matched-pair group comparison workflow with permutation tests and BH-FDR.
- `compareCrossFrequencyPhaseAmplitudeGrangerProtocolDifferenceBK1.m`: Compares two protocols directly using only the matched pairs with complete data in both protocols.
- `runAnalysis2CrossFrequencyPhaseAmplitudeGrangerBK1.m`: Produces compact secondary MAT summaries for protocol consistency, effect sizes, heterogeneity, and raw-vs-normalized inference.
- `runNetFlowAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1.m`: Recomputes ROI-wise and band-wise net flow from saved CF-GC results and saves dedicated overview figures.

### Helper Package: `+cfGCUtils`

- `addFigureTitle.m`: Adds a figure-wide title while staying compatible with MATLAB versions that do not support `sgtitle`.
- `bhFDRByPair.m`: Applies BH-FDR separately within each phase/amplitude band pair.
- `blueWhiteRed.m`: Returns the diverging colormap used for signed difference plots.
- `computeNetFlow.m`: Converts a directed matrix into ROI net-flow values.
- `createBandPairGridFigure.m`: Creates one tiled heatmap figure covering all selected phase-amplitude band pairs.
- `createFigure.m`: Creates hidden white figures with standard sizing defaults.
- `ensureFolder.m`: Creates output folders when needed.
- `getBK1SubjectGroups.m`: Loads BK1 subjects and labels them as meditator, control, or unknown.
- `getDefaultProtocolNameList.m`: Returns the default BK1 protocol order.
- `getMatchedSubjectPairsBK1.m`: Loads the BK1 matched-pair table.
- `loadConditionResults.m`: Loads one saved CF-GC MAT file for a given subject/protocol/condition.
- `makeAnalysisTag.m`: Builds stable folder tags from human-readable labels plus hashed settings.
- `makeOptionsHash.m`: Generates the settings hash used for cache-safe output folders.
- `maxFiniteValue.m`: Returns the largest finite array value.
- `meanAcrossPairs.m`: Collapses a band-pair stack by averaging across phase/amplitude combinations.
- `minAcrossPairs.m`: Collapses a band-pair stack by taking the minimum across combinations.
- `minFiniteValue.m`: Returns the smallest finite array value.
- `normalizeList.m`: Normalizes string, char, or cell inputs into a standard cell-array format.
- `normalizeSubjectGrid.m`: Z-scores a subject matrix using only finite entries.
- `pairedLabelSwapMatchedGroups.m`: Runs matched-pair permutation tests for meditator-versus-control contrasts.
- `pairedLabelSwapPrePost.m`: Runs paired permutation tests for PRE-versus-POST contrasts.
- `plotHeatmap.m`: Draws consistently formatted ROI heatmaps.
- `resolveColormap.m`: Accepts a colormap name, function handle, or raw matrix and resolves it to a usable colormap.
- `resolveCommonCLim.m`: Computes shared color limits across one or more matrices.
- `resolveDiffCLim.m`: Computes symmetric color limits for signed difference maps.
- `resolvePairSelectionBK1.m`: Converts user pair input into a validated BK1 matched-pair list.
- `resolveProtocolNameList.m`: Expands `'all'` and validates protocol-name inputs.
- `resolveSubjectSelection.m`: Filters subjects to those with complete CF-GC files for the requested protocols.
- `runMultiplePairComparisonProtocol.m`: Executes the per-protocol worker used by the multi-pair comparison entry point.
- `runSinglePairComparisonProtocol.m`: Executes the per-protocol worker used by the single-pair comparison entry point.
- `runSubjectVisualizationProtocol.m`: Executes the per-protocol worker used by the single-subject visualization entry point.
- `sanitizePathToken.m`: Converts free text into safe file/folder tokens.
- `saveAnalysisSummary.m`: Saves an `analysisResult` structure after creating its parent folder.
- `saveFigureQuietly.m`: Saves and closes figures without leaving GUI clutter behind.
- `sumAcrossPairs.m`: Collapses a band-pair stack by summing across combinations.
- `tryLoadCachedAnalysis.m`: Reloads cached analysis outputs when the expected files are still present.
- `unpackBandPairGrid.m`: Extracts the selected CF-GC subgrid and band metadata from a saved result.
- `unpairedLabelShuffle.m`: Runs a permutation test for independent-group contrasts.
- `validateResultCompatibility.m`: Confirms that two CF-GC result structs can be compared safely.


