function analysisResult = compareCrossFrequencyPhaseAmplitudeGrangerMultiplePairsBK1(pairSelection,protocolNameList,badEyeCondition,badTrialVersion,cfgIn)
% compareCrossFrequencyPhaseAmplitudeGrangerMultiplePairsBK1
%
% Compare matched meditator-control pairs across one or more protocols.
% This entry point loads saved PRE/POST CF-GC results, aggregates the
% selected BK1 pairs, runs permutation-based inference, and saves summary
% figures plus MAT files under analysis/comparison/multiple_pairs.
%
% pairSelection can be:
%   'all'
%   a numeric vector of matched-pair indices
%   a cell array of explicit pairs such as
%   {{'013AR','064PK'}; {'019CKa','022SSP'}}

if ~exist('pairSelection','var') || isempty(pairSelection);     pairSelection = 'all'; end
if ~exist('protocolNameList','var') || isempty(protocolNameList); protocolNameList = 'all'; end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('cfgIn','var') || isempty(cfgIn);                     cfgIn = struct; end

cfg = applyDefaultCfg(cfgIn);
protocolNameList = cfGCUtils.resolveProtocolNameList(protocolNameList);
projectRoot = fileparts(fileparts(mfilename('fullpath')));
[pairList,selectionLabel] = cfGCUtils.resolvePairSelectionBK1(pairSelection,projectRoot);

optionsText = sprintf('multiple_pairs|pairs=%s|protocols=%s|eye=%s|trial=%s|perm=%d|minobs=%d|minpairs=%d|fdr=%.4f|overview=%d|save=%d', ...
    strjoin({pairList.pairLabel},','),strjoin(protocolNameList,','),badEyeCondition,badTrialVersion, ...
    cfg.numPermutations,cfg.minObservationCount,cfg.minValidSubjectsPerCell,cfg.fdrAlpha,cfg.makeOverviewFigure,cfg.saveOutputFlag);
analysisTag = cfGCUtils.makeAnalysisTag(selectionLabel,optionsText);
analysisRoot = buildAnalysisRoot(cfg.analysisFolder,cfg.analysisSubfolder,analysisTag);
summaryFileName = fullfile(analysisRoot,'multiple_pairs_comparison_summary.mat');

[cachedResult,wasLoaded] = cfGCUtils.tryLoadCachedAnalysis(summaryFileName,cfg.saveOutputFlag && ~cfg.forceRebuild);
if wasLoaded && ~cfg.forceRebuild
    analysisResult = cachedResult;
    fprintf('Loaded cached multi-pair comparison:\n%s\n',summaryFileName);
    return;
end

analysisResult = struct;
analysisResult.mainOption = 'comparison';
analysisResult.subOption = 'multiple_pairs';
analysisResult.pairList = pairList;
analysisResult.protocolNameList = protocolNameList;
analysisResult.badEyeCondition = badEyeCondition;
analysisResult.badTrialVersion = badTrialVersion;
analysisResult.analysisFolder = analysisRoot;
analysisResult.summaryFileName = summaryFileName;
analysisResult.protocolResults = repmat(struct('protocolName','','result',struct(),'summaryFileName',''),numel(protocolNameList),1);
analysisResult.outputManifest = {};

for iProtocol = 1:numel(protocolNameList)
    protocolName = protocolNameList{iProtocol};
    protocolRoot = fullfile(analysisRoot,protocolName);
    protocolCfg = cfg;
    protocolCfg.figureFolder = fullfile(protocolRoot,'figures');

    protocolResult = cfGCUtils.runMultiplePairComparisonProtocol(pairList,protocolName,badEyeCondition,badTrialVersion,protocolCfg);
    protocolSummaryFileName = fullfile(protocolRoot,[protocolName '_multiple_pairs_comparison.mat']);
    cfGCUtils.saveAnalysisSummary(protocolSummaryFileName,protocolResult);

    analysisResult.protocolResults(iProtocol).protocolName = protocolName;
    analysisResult.protocolResults(iProtocol).result = protocolResult;
    analysisResult.protocolResults(iProtocol).summaryFileName = protocolSummaryFileName;
    analysisResult.outputManifest = [analysisResult.outputManifest; protocolResult.outputManifest; {protocolSummaryFileName}]; %#ok<AGROW>
end

cfGCUtils.saveAnalysisSummary(summaryFileName,analysisResult);
analysisResult.outputManifest{end+1,1} = summaryFileName;
cfGCUtils.saveAnalysisSummary(summaryFileName,analysisResult);
end

function cfg = applyDefaultCfg(cfgIn)
thisFolder = fileparts(mfilename('fullpath'));

cfg = cfgIn;
cfg = setDefault(cfg,'dataFolder',fullfile(thisFolder,'savedDataCrossFreqGranger'));
cfg = setDefault(cfg,'analysisFolder',fullfile(thisFolder,'analysis'));
cfg = setDefault(cfg,'analysisSubfolder',fullfile('comparison','multiple_pairs'));
cfg = setDefault(cfg,'numPermutations',5000);
cfg = setDefault(cfg,'minObservationCount',20);
cfg = setDefault(cfg,'minValidSubjectsPerCell',3);
cfg = setDefault(cfg,'fdrAlpha',0.05);
cfg = setDefault(cfg,'commonColorLimits',[]);
cfg = setDefault(cfg,'diffColorLimits',[]);
cfg = setDefault(cfg,'makeOverviewFigure',1);
cfg = setDefault(cfg,'saveOutputFlag',1);
cfg = setDefault(cfg,'forceRebuild',0);
end

function analysisRoot = buildAnalysisRoot(baseFolder,subFolder,analysisTag)
if isempty(subFolder)
    analysisRoot = fullfile(baseFolder,analysisTag);
else
    analysisRoot = fullfile(baseFolder,subFolder,analysisTag);
end
end

function s = setDefault(s,fieldName,defaultValue)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end
