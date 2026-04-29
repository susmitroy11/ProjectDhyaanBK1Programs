function analysisResult = visualizeCrossFrequencyPhaseAmplitudeGrangerSingleSubjectBK1(subjectName,protocolNameList,badEyeCondition,badTrialVersion,cfgIn)
% visualizeCrossFrequencyPhaseAmplitudeGrangerSingleSubjectBK1
%
% Descriptively visualize CF-GC results for one subject.
% This entry point loads saved PRE/POST CF-GC results for a single BK1
% subject, creates band-pair heatmaps plus overview figures, and saves the
% outputs under analysis/visualization/single_subject.

if ~exist('subjectName','var') || isempty(subjectName);         subjectName = '019CKa'; end
if ~exist('protocolNameList','var') || isempty(protocolNameList); protocolNameList = 'all'; end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('cfgIn','var') || isempty(cfgIn);                     cfgIn = struct; end

cfg = applyDefaultCfg(cfgIn);
protocolNameList = cfGCUtils.resolveProtocolNameList(protocolNameList);

optionsText = sprintf('single_subject|subject=%s|protocols=%s|eye=%s|trial=%s|sample=%d|overview=%d|samplefig=%d|save=%d', ...
    subjectName,strjoin(protocolNameList,','),badEyeCondition,badTrialVersion, ...
    cfg.sampleCountThreshold,cfg.showOverviewFigure,cfg.showSampleCountFigure,cfg.saveOutputFlag);
analysisTag = cfGCUtils.makeAnalysisTag(subjectName,optionsText);
analysisRoot = buildAnalysisRoot(cfg.analysisFolder,cfg.analysisSubfolder,analysisTag);
summaryFileName = fullfile(analysisRoot,'single_subject_visualization_summary.mat');

[cachedResult,wasLoaded] = cfGCUtils.tryLoadCachedAnalysis(summaryFileName,cfg.saveOutputFlag && ~cfg.forceRebuild);
if wasLoaded && ~cfg.forceRebuild
    analysisResult = cachedResult;
    fprintf('Loaded cached single-subject visualization:\n%s\n',summaryFileName);
    return;
end

analysisResult = struct;
analysisResult.mainOption = 'visualization';
analysisResult.subOption = 'single_subject';
analysisResult.subjectName = subjectName;
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

    protocolResult = cfGCUtils.runSubjectVisualizationProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,protocolCfg);
    protocolSummaryFileName = fullfile(protocolRoot,[subjectName '_' protocolName '_single_subject_visualization.mat']);
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
cfg = setDefault(cfg,'analysisSubfolder',fullfile('visualization','single_subject'));
cfg = setDefault(cfg,'sampleCountThreshold',20);
cfg = setDefault(cfg,'commonColorLimits',[]);
cfg = setDefault(cfg,'diffColorLimits',[]);
cfg = setDefault(cfg,'showSampleCountFigure',1);
cfg = setDefault(cfg,'showOverviewFigure',1);
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
