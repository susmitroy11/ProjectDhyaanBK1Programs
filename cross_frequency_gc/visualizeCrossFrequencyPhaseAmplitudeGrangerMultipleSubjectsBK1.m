function analysisResult = visualizeCrossFrequencyPhaseAmplitudeGrangerMultipleSubjectsBK1(subjectSelection,protocolNameList,badEyeCondition,badTrialVersion,cfgIn)
% visualizeCrossFrequencyPhaseAmplitudeGrangerMultipleSubjectsBK1
%
% Descriptively visualize CF-GC results for multiple subjects.
% This entry point loops over a validated subject list, calls the
% single-subject visualization worker for each subject, and saves one
% combined summary under analysis/visualization/multiple_subjects.
%
% subjectSelection can be 'all' or an explicit list such as
% {'013AR','064PK'}.

if ~exist('subjectSelection','var') || isempty(subjectSelection); subjectSelection = 'all'; end
if ~exist('protocolNameList','var') || isempty(protocolNameList); protocolNameList = 'all'; end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition);   badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion);   badTrialVersion = 'v8'; end
if ~exist('cfgIn','var') || isempty(cfgIn);                       cfgIn = struct; end

cfg = applyDefaultCfg(cfgIn);
protocolNameList = cfGCUtils.resolveProtocolNameList(protocolNameList);
[subjectNameList,selectionLabel] = cfGCUtils.resolveSubjectSelection(subjectSelection,cfg.dataFolder,protocolNameList,badEyeCondition,badTrialVersion);

optionsText = sprintf('multiple_subjects|subjects=%s|protocols=%s|eye=%s|trial=%s|sample=%d|overview=%d|samplefig=%d|save=%d', ...
    strjoin(subjectNameList,','),strjoin(protocolNameList,','),badEyeCondition,badTrialVersion, ...
    cfg.sampleCountThreshold,cfg.showOverviewFigure,cfg.showSampleCountFigure,cfg.saveOutputFlag);
analysisTag = cfGCUtils.makeAnalysisTag(selectionLabel,optionsText);
analysisRoot = fullfile(cfg.analysisFolder,'visualization','multiple_subjects',analysisTag);
summaryFileName = fullfile(analysisRoot,'multiple_subjects_visualization_summary.mat');

[cachedResult,wasLoaded] = cfGCUtils.tryLoadCachedAnalysis(summaryFileName,cfg.saveOutputFlag && ~cfg.forceRebuild);
if wasLoaded && ~cfg.forceRebuild
    analysisResult = cachedResult;
    fprintf('Loaded cached multi-subject visualization:\n%s\n',summaryFileName);
    return;
end

analysisResult = struct;
analysisResult.mainOption = 'visualization';
analysisResult.subOption = 'multiple_subjects';
analysisResult.subjectNameList = subjectNameList;
analysisResult.protocolNameList = protocolNameList;
analysisResult.badEyeCondition = badEyeCondition;
analysisResult.badTrialVersion = badTrialVersion;
analysisResult.analysisFolder = analysisRoot;
analysisResult.summaryFileName = summaryFileName;
analysisResult.subjectResults = repmat(struct('subjectName','','protocolResults',struct([]),'summaryFileName',''),numel(subjectNameList),1);
analysisResult.outputManifest = {};

for iSubject = 1:numel(subjectNameList)
    subjectName = subjectNameList{iSubject};
    subjectCfg = cfg;
    subjectCfg.analysisFolder = fullfile(analysisRoot,subjectName);
    subjectCfg.analysisSubfolder = '';
    subjectResult = visualizeCrossFrequencyPhaseAmplitudeGrangerSingleSubjectBK1(subjectName,protocolNameList,badEyeCondition,badTrialVersion,subjectCfg);
    subjectSummaryFileName = fullfile(subjectCfg.analysisFolder,'subject_summary.mat');
    cfGCUtils.saveAnalysisSummary(subjectSummaryFileName,subjectResult);

    analysisResult.subjectResults(iSubject).subjectName = subjectName;
    analysisResult.subjectResults(iSubject).protocolResults = subjectResult.protocolResults;
    analysisResult.subjectResults(iSubject).summaryFileName = subjectSummaryFileName;
    analysisResult.outputManifest = [analysisResult.outputManifest; subjectResult.outputManifest; {subjectSummaryFileName}]; %#ok<AGROW>
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
cfg = setDefault(cfg,'analysisSubfolder',fullfile('visualization','multiple_subjects'));
cfg = setDefault(cfg,'sampleCountThreshold',20);
cfg = setDefault(cfg,'commonColorLimits',[]);
cfg = setDefault(cfg,'diffColorLimits',[]);
cfg = setDefault(cfg,'showSampleCountFigure',1);
cfg = setDefault(cfg,'showOverviewFigure',1);
cfg = setDefault(cfg,'saveOutputFlag',1);
cfg = setDefault(cfg,'forceRebuild',0);
end

function s = setDefault(s,fieldName,defaultValue)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end
