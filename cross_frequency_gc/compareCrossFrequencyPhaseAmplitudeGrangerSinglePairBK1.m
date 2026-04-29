function analysisResult = compareCrossFrequencyPhaseAmplitudeGrangerSinglePairBK1(pairInput,protocolNameList,badEyeCondition,badTrialVersion,cfgIn)
% compareCrossFrequencyPhaseAmplitudeGrangerSinglePairBK1
%
% Descriptively compare one matched meditator-control pair.
% This entry point loads saved PRE/POST CF-GC results for one BK1 matched
% pair, generates figure panels for each requested protocol, and saves a
% compact summary structure under analysis/comparison/single_pair.
%
% pairInput can be a numeric pair index or an explicit pair such as
% {'013AR','064PK'}.

if ~exist('pairInput','var') || isempty(pairInput);             pairInput = 1; end
if ~exist('protocolNameList','var') || isempty(protocolNameList); protocolNameList = 'all'; end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('cfgIn','var') || isempty(cfgIn);                     cfgIn = struct; end

cfg = applyDefaultCfg(cfgIn);
protocolNameList = cfGCUtils.resolveProtocolNameList(protocolNameList);
projectRoot = fileparts(fileparts(mfilename('fullpath')));
[pairList,~] = cfGCUtils.resolvePairSelectionBK1(pairInput,projectRoot);
if numel(pairList) ~= 1
    error('Single-pair comparison expects exactly one pair.');
end
pairSelection = pairList(1);

optionsText = sprintf('single_pair|pair=%s|protocols=%s|eye=%s|trial=%s|minobs=%d|overview=%d|save=%d', ...
    pairSelection.pairLabel,strjoin(protocolNameList,','),badEyeCondition,badTrialVersion, ...
    cfg.minObservationCount,cfg.makeOverviewFigure,cfg.saveOutputFlag);
analysisTag = cfGCUtils.makeAnalysisTag(pairSelection.pairLabel,optionsText);
analysisRoot = buildAnalysisRoot(cfg.analysisFolder,cfg.analysisSubfolder,analysisTag);
summaryFileName = fullfile(analysisRoot,'single_pair_comparison_summary.mat');

[cachedResult,wasLoaded] = cfGCUtils.tryLoadCachedAnalysis(summaryFileName,cfg.saveOutputFlag && ~cfg.forceRebuild);
if wasLoaded && ~cfg.forceRebuild
    analysisResult = cachedResult;
    fprintf('Loaded cached single-pair comparison:\n%s\n',summaryFileName);
    return;
end

analysisResult = struct;
analysisResult.mainOption = 'comparison';
analysisResult.subOption = 'single_pair';
analysisResult.pairSelection = pairSelection;
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

    protocolResult = cfGCUtils.runSinglePairComparisonProtocol(pairSelection,protocolName,badEyeCondition,badTrialVersion,protocolCfg);
    protocolSummaryFileName = fullfile(protocolRoot,[protocolName '_' pairSelection.pairLabel '_single_pair_comparison.mat']);
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
cfg = setDefault(cfg,'analysisSubfolder',fullfile('comparison','single_pair'));
cfg = setDefault(cfg,'minObservationCount',20);
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
