function saveGrangerData(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,stRange,gcParams)
% saveGrangerData Save multivariate spectral Granger causality for BK1 EEG.
%
% Usage:
%   saveGrangerData(subjectName,protocolNameList,badEyeCondition,...
%       badTrialVersion,ftDataFolder,stRange,gcParams)
%
% Inputs:
%   subjectName       : BK1 subject code, e.g. '028HB'
%   protocolNameList  : cell array of protocol names, e.g. {'M1' 'M2'}
%   badEyeCondition   : 'ep' or 'wo'
%   badTrialVersion   : bad trial version string, e.g. 'v8'
%   ftDataFolder      : folder containing per-subject FieldTrip data
%   stRange           : stimulus/post window in seconds, e.g. [0.25 1.25]
%   gcParams          : optional struct with fields:
%       .modelOrder       default 10
%       .mvarToolbox      default 'biosig'
%       .biosigMvarMethod default 2
%       .maxFreq          default 100
%       .freqStep         default 1
%       .feedback         default 'none'
%       .electrodeList    default [] (use all non-bad electrodes)
%       .flatThreshold    default 1e-6
%       .minSampleRatio   default 3
%       .savedDataFolder  default <this folder>/savedDataGranger
%
% Output:
%   Saves one file per protocol as:
%   <savedDataFolder>/<subject>/<protocol>_<badEyeCondition>_<badTrialVersion>_granger.mat

if ~exist('stRange','var') || isempty(stRange)
    stRange = [0.25 1.25];
end
if ~exist('gcParams','var') || isempty(gcParams)
    gcParams = struct;
end

gcParams = applyDefaultParams(gcParams);
ensureLocalFieldTripOnPath();

for iProtocol = 1:numel(protocolNameList)
    protocolName = protocolNameList{iProtocol};
    saveGrangerDataSingleProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,ftDataFolder,stRange,gcParams);
end
end

function saveGrangerDataSingleProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,ftDataFolder,stRange,gcParams)

savedDataFolder = fullfile(gcParams.savedDataFolder,subjectName);
ensureFolder(savedDataFolder);

analysisFileName = fullfile(savedDataFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_granger.mat']);
ftDataFileName = fullfile(ftDataFolder,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '.mat']);

numGoodTrials = 0;
label = {};
grangerPre = [];
grangerPost = [];
freqPre = [];
freqPost = [];
metaPre = struct;
metaPost = struct;
gcParamsUsed = gcParams;

if ~exist(ftDataFileName,'file')
    metaPre.status = 'missing_ft_data';
    metaPost.status = 'missing_ft_data';
    save(analysisFileName,'grangerPre','grangerPost','freqPre','freqPost','numGoodTrials','label','metaPre','metaPost','gcParamsUsed','stRange');
    return;
end

tmpData = load(ftDataFileName);
numGoodTrials = tmpData.numGoodTrials;

if numGoodTrials <= 0 || isempty(tmpData.data)
    if isfield(tmpData,'data') && isfield(tmpData.data,'label')
        label = tmpData.data.label;
    end
    metaPre.status = 'no_good_trials';
    metaPost.status = 'no_good_trials';
    save(analysisFileName,'grangerPre','grangerPost','freqPre','freqPost','numGoodTrials','label','metaPre','metaPost','gcParamsUsed','stRange');
    return;
end

data = tmpData.data;
label = data.label;

cfg = [];
cfg.toilim = [(-diff(stRange) + 1/data.fsample) 0];
dataPre = ft_redefinetrial(cfg,data);

cfg = [];
cfg.toilim = [(stRange(1) + 1/data.fsample) stRange(2)];
dataPost = ft_redefinetrial(cfg,data);

[grangerPre,freqPre,metaPre] = getGrangerThisCondition(dataPre,gcParams);
[grangerPost,freqPost,metaPost] = getGrangerThisCondition(dataPost,gcParams);

save(analysisFileName,'grangerPre','grangerPost','freqPre','freqPost','numGoodTrials','label','metaPre','metaPost','gcParamsUsed','stRange');
end

function [grangerFull,freqVals,meta] = getGrangerThisCondition(data,gcParams)

numAllChannels = numel(data.label);
freqVals = [];
grangerFull = [];

meta = struct;
meta.status = 'not_run';
meta.modelOrder = gcParams.modelOrder;
meta.mvarToolbox = gcParams.mvarToolbox;
meta.maxFreq = gcParams.maxFreq;
meta.freqStep = gcParams.freqStep;
meta.badElectrodes = [];
meta.flatElectrodes = [];
meta.selectedElectrodes = [];
meta.selectedLabels = {};
meta.numChannels = 0;
meta.numTrials = numel(data.trial);
meta.samplesPerParameter = NaN;
meta.minSampleRatio = gcParams.minSampleRatio;
meta.lowSampleRatioFlag = false;

badElecs = [];
if isfield(data,'badElecs') && ~isempty(data.badElecs)
    badElecs = unique(data.badElecs(:)');
    badElecs = badElecs(badElecs >= 1 & badElecs <= numAllChannels);
end
meta.badElectrodes = badElecs;

selectedIndices = setdiff(1:numAllChannels,badElecs,'stable');
if ~isempty(gcParams.electrodeList)
    requestedElecs = unique(gcParams.electrodeList(:)');
    requestedElecs = requestedElecs(requestedElecs >= 1 & requestedElecs <= numAllChannels);
    selectedIndices = intersect(selectedIndices,requestedElecs,'stable');
end

flatElecs = findFlatChannels(data,selectedIndices,gcParams.flatThreshold);
selectedIndices = setdiff(selectedIndices,flatElecs,'stable');
meta.flatElectrodes = flatElecs;

if numel(selectedIndices) < 2
    meta.status = 'too_few_channels';
    return;
end

cfg = [];
cfg.channel = data.label(selectedIndices);
dataUse = ft_selectdata(cfg,data);

totalSamples = sum(cellfun(@(x) size(x,2),dataUse.trial));
numChannels = numel(dataUse.label);
numParameters = max(1,(numChannels^2) * gcParams.modelOrder);
meta.samplesPerParameter = totalSamples / numParameters;
meta.lowSampleRatioFlag = meta.samplesPerParameter < gcParams.minSampleRatio;

meta.selectedElectrodes = selectedIndices;
meta.selectedLabels = data.label(selectedIndices);
meta.numChannels = numChannels;

try
    cfg = [];
    cfg.method = gcParams.mvarToolbox;
    cfg.order = gcParams.modelOrder;
    cfg.keeptrials = 'no';
    cfg.demean = 'yes';
    cfg.feedback = gcParams.feedback;
    if strcmpi(gcParams.mvarToolbox,'biosig')
        cfg.mvarmethod = gcParams.biosigMvarMethod;
    end
    mvarData = ft_mvaranalysis(cfg,dataUse);

    cfg = [];
    cfg.method = 'mvar';
    cfg.foi = 0:gcParams.freqStep:min(gcParams.maxFreq,floor(dataUse.fsample/2));
    cfg.feedback = gcParams.feedback;
    freqData = ft_freqanalysis(cfg,mvarData);

    cfg = [];
    cfg.method = 'granger';
    cfg.feedback = gcParams.feedback;
    grangerStat = ft_connectivityanalysis(cfg,freqData);

    freqVals = grangerStat.freq;
    grangerFull = nan(numAllChannels,numAllChannels,numel(freqVals));
    grangerFull(selectedIndices,selectedIndices,:) = grangerStat.grangerspctrm;

    meta.status = 'ok';
catch ME
    meta.status = 'error';
    meta.errorMessage = ME.message;
end
end

function flatElecs = findFlatChannels(data,selectedIndices,flatThreshold)

flatElecs = [];
if isempty(selectedIndices) || isempty(data.trial)
    return;
end

allSamples = cat(2,data.trial{:});
channelStd = std(allSamples(selectedIndices,:),0,2);
flatElecs = selectedIndices(channelStd <= flatThreshold);
end

function gcParams = applyDefaultParams(gcParams)

    defaultSaveFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataGranger');

gcParams = setDefault(gcParams,'modelOrder',10);
gcParams = setDefault(gcParams,'mvarToolbox','biosig');
gcParams = setDefault(gcParams,'biosigMvarMethod',2);
gcParams = setDefault(gcParams,'maxFreq',100);
gcParams = setDefault(gcParams,'freqStep',1);
gcParams = setDefault(gcParams,'feedback','none');
gcParams = setDefault(gcParams,'electrodeList',[]);
gcParams = setDefault(gcParams,'flatThreshold',1e-6);
gcParams = setDefault(gcParams,'minSampleRatio',3);
gcParams = setDefault(gcParams,'savedDataFolder',defaultSaveFolder);
end

function s = setDefault(s,fieldName,defaultValue)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end

function ensureLocalFieldTripOnPath()
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
fieldTripPath = fullfile(projectRoot,'fieldtrip-20260211');
addpath(fieldTripPath);
ft_defaults;
end
