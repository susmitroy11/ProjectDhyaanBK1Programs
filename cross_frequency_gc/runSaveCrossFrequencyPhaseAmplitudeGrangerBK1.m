function batchResults = runSaveCrossFrequencyPhaseAmplitudeGrangerBK1(protocolNameList,conditionTypeList,badEyeCondition,badTrialVersion,cfgIn)
% runSaveCrossFrequencyPhaseAmplitudeGrangerBK1
%
% Batch runner for subject-level CF-GC computation.
% Use this function to generate the savedDataCrossFreqGranger MAT files
% that the visualization and comparison scripts expect.
%
% This calls:
%   computeCrossFrequencyPhaseAmplitudeGrangerBK1
%
% for all good BK1 subjects, all requested protocols, and all requested
% condition windows, and saves one MAT file per subject/protocol/condition.
%
% The underlying computation is restricted to the same forward-only band
% pairs as the reference script:
%   phase 4-8, 8-12, 12-16 Hz  -> amplitude 32-36, 36-40, 40-44 Hz
% across all 6x6 ROI combinations.
%
% Examples:
%   batchResults = runSaveCrossFrequencyPhaseAmplitudeGrangerBK1;
%
%   protocolNameList = {'M1','M2'};
%   conditionTypeList = {'post'};
%   cfg = struct;
%   cfg.modelOrder = 10;
%   batchResults = runSaveCrossFrequencyPhaseAmplitudeGrangerBK1(protocolNameList,conditionTypeList,'ep','v8',cfg);

if ~exist('protocolNameList','var') || isempty(protocolNameList)
    protocolNameList = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
end
if ~exist('conditionTypeList','var') || isempty(conditionTypeList)
    conditionTypeList = {'pre','post'};
end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition)
    badEyeCondition = 'ep';
end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion)
    badTrialVersion = 'v8';
end
if ~exist('cfgIn','var') || isempty(cfgIn)
    cfgIn = struct;
end
cfgIn = applyRunnerDefaults(cfgIn);

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(projectRoot,'ProjectDhyaanBK1Programs')));
addpath(genpath(fullfile(projectRoot,'CommonPrograms')));
addpath(genpath(fullfile(projectRoot,'Montages')));
addpath(fullfile(projectRoot,'fieldtrip-20260211'));
ft_defaults;

infoFolder = fullfile(projectRoot,'ProjectDhyaanBK1Programs','commonAnalysisCodes','informationFiles');
addpath(infoFolder);
cleanupObj = onCleanup(@() rmpath(infoFolder)); %#ok<NASGU>

[goodSubjectList,~,~] = getGoodSubjectsBK1();
goodSubjectList = normalizeList(goodSubjectList);

numSubjects = numel(goodSubjectList);
numProtocols = numel(protocolNameList);
numConditions = numel(conditionTypeList);
numJobs = numSubjects * numProtocols * numConditions;
parallelAvailable = false;

if cfgIn.useParallel || cfgIn.parallelizeJobs
    parallelAvailable = ensureParallelPool(cfgIn.numWorkers);
end

batchResults = struct;
batchResults.projectRoot = projectRoot;
batchResults.protocolNameList = protocolNameList;
batchResults.conditionTypeList = conditionTypeList;
batchResults.badEyeCondition = badEyeCondition;
batchResults.badTrialVersion = badTrialVersion;
batchResults.cfgIn = cfgIn;
batchResults.parallelAvailable = parallelAvailable;
batchResults.goodSubjectList = goodSubjectList;
batchResults.jobTable = repmat(struct( ...
    'subjectName','', ...
    'protocolName','', ...
    'conditionType','', ...
    'status','pending', ...
    'outputFileName','', ...
    'errorMessage',''), numJobs, 1);

jobSpecs = batchResults.jobTable;
jobCounter = 0;
for iSubject = 1:numSubjects
    subjectName = goodSubjectList{iSubject};
    for iProtocol = 1:numProtocols
        protocolName = protocolNameList{iProtocol};
        for iCondition = 1:numConditions
            conditionType = conditionTypeList{iCondition};
            jobCounter = jobCounter + 1;
            jobSpecs(jobCounter).subjectName = subjectName;
            jobSpecs(jobCounter).protocolName = protocolName;
            jobSpecs(jobCounter).conditionType = conditionType;
        end
    end
end

savedDataRoot = getSavedDataRoot(cfgIn, fileparts(mfilename('fullpath')));
analysisNotes = '';

if cfgIn.parallelizeJobs && parallelAvailable
    workerCfg = cfgIn;
    workerCfg.useParallel = 0;
    parfor iJob = 1:numJobs
        jobSpecs(iJob) = runOneJob(jobSpecs(iJob), badEyeCondition, badTrialVersion, workerCfg, savedDataRoot);
    end
else
    for iJob = 1:numJobs
        jobSpecs(iJob) = runOneJob(jobSpecs(iJob), badEyeCondition, badTrialVersion, cfgIn, savedDataRoot);
    end
end

batchResults.jobTable = jobSpecs;
for iJob = 1:numJobs
    if strcmp(jobSpecs(iJob).status,'ok')
        analysisNotes = ['Exploratory regional phase-amplitude Granger batch completed ', ...
            '(forward-only, restricted low-phase/high-amplitude band pairs).'];
        break;
    end
end
if ~isempty(analysisNotes)
    batchResults.analysisNotes = analysisNotes;
end

statusList = {batchResults.jobTable.status};
batchResults.numJobs = numJobs;
batchResults.numSuccess = sum(strcmp(statusList,'ok'));
batchResults.numError = sum(strcmp(statusList,'error'));

summaryFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataCrossFreqGranger');
ensureFolder(summaryFolder);
summaryFileName = fullfile(summaryFolder,['batchSummary_' badEyeCondition '_' badTrialVersion '.mat']);
save(summaryFileName,'batchResults','-v7.3');
batchResults.summaryFileName = summaryFileName;

fprintf('\nCompleted batch run: %d/%d successful, %d failed.\n', ...
    batchResults.numSuccess,batchResults.numJobs,batchResults.numError);
fprintf('Batch summary saved to:\n%s\n',summaryFileName);
end

function cfgOut = applyRunnerDefaults(cfgIn)
cfgOut = cfgIn;
cfgOut = setDefault(cfgOut,'useParallel',1);
cfgOut = setDefault(cfgOut,'parallelizeJobs',0);
cfgOut = setDefault(cfgOut,'numWorkers',[]);
end

function s = setDefault(s,fieldName,defaultValue)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end

function tf = ensureParallelPool(numWorkers)
tf = false;
if exist('parpool','file') ~= 2 || exist('gcp','file') ~= 2
    return;
end
try
    poolObj = gcp('nocreate');
    if isempty(poolObj)
        if isempty(numWorkers)
            parpool('local');
        else
            parpool('local',numWorkers);
        end
    end
    tf = true;
catch ME
    fprintf(2,'Could not start parallel pool: %s\n',ME.message);
    tf = false;
end
end

function savedDataRoot = getSavedDataRoot(cfgIn,codeFolder)
if isfield(cfgIn,'savedDataFolder') && ~isempty(cfgIn.savedDataFolder)
    savedDataRoot = cfgIn.savedDataFolder;
else
    savedDataRoot = fullfile(codeFolder,'savedDataCrossFreqGranger');
end
end

function jobResult = runOneJob(jobResult,badEyeCondition,badTrialVersion,cfgIn,savedDataRoot)
subjectName = jobResult.subjectName;
protocolName = jobResult.protocolName;
conditionType = jobResult.conditionType;

fprintf('Running cross-frequency phase-amplitude GC | subject %s | protocol %s | %s\n', ...
    subjectName, protocolName, upper(conditionType));

try
    computeCrossFrequencyPhaseAmplitudeGrangerBK1( ...
        subjectName,protocolName,conditionType,badEyeCondition,badTrialVersion,cfgIn);

    jobResult.status = 'ok';
    jobResult.outputFileName = fullfile(savedDataRoot,subjectName, ...
        [protocolName '_' conditionType '_' badEyeCondition '_' badTrialVersion '_crossFreqPhaseAmpGC.mat']);
catch ME
    jobResult.status = 'error';
    jobResult.errorMessage = ME.message;
    fprintf(2,'Error for %s %s %s: %s\n',subjectName,protocolName,conditionType,ME.message);
end
end

function out = normalizeList(in)
if iscell(in)
    out = in(:);
elseif isstring(in)
    out = cellstr(in(:));
else
    out = cellstr(in);
end
end

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
