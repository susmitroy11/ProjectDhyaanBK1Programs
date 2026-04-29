# Spectral Granger Causality

## Overview

This folder contains MATLAB utilities for plotting BK1 spectral Granger flow as scalp topographies. The main implementation is `plotGrangerOutflowTopoplotBK1.m`, and the other two files are convenience wrappers for inflow and net-flow views.

The code expects saved spectral Granger files in `savedDataGranger/` and uses BK1 subject-group metadata plus the actiCap montage and EEGLAB `topoplot`.

## Typical Usage

Examples:

```matlab
plotGrangerOutflowTopoplotBK1('028HB','M1',[30 80],'diff');
plotGrangerInflowTopoplotBK1('Meditator','M1',[30 80],'diff','ep','v8',1);
plotGrangerNetFlowTopoplotBK1('BothGroups','all',[30 80],'all','ep','v8',1);
```

Useful input conventions:

- `selection`: a subject code such as `'028HB'`, or a group like `'Meditator'`, `'Control'`, or `'All'`
- `protocolName`: a specific protocol such as `'M1'`, or `'all'`
- `conditionType`: `'pre'`, `'post'`, `'diff'`, or `'all'`
- `freqRange`: two-element frequency window in Hz

## File Guide

- `plotGrangerOutflowTopoplotBK1.m`: Main plotting function for outflow, inflow, or net flow across one subject or grouped subjects.
- `plotGrangerInflowTopoplotBK1.m`: Wrapper that calls the main function with `measureType = 'inflow'`.
- `plotGrangerNetFlowTopoplotBK1.m`: Wrapper that calls the main function with `measureType = 'netflow'`.

## Requirements

- MATLAB with access to the BK1 project folders
- Saved spectral Granger MAT files under `savedDataGranger/`
- Montage files for `actiCap64_UOL`
- EEGLAB `topoplot` available directly or through the bundled FieldTrip external EEGLAB copy

