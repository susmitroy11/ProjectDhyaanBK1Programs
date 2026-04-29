function results = computeCrossFrequencyPhaseAmplitudeGrangerBK1(subjectName,protocolName,conditionType,badEyeCondition,badTrialVersion,cfgIn)
% computeCrossFrequencyPhaseAmplitudeGrangerBK1
%
% Compute subject-level cross-frequency phase-amplitude Granger matrices.
% This is the core BK1 CF-GC worker used by the batch runner and all later
% visualization/comparison scripts.
%
% Important methodological note:
%   Doing MVAR/Granger directly on wrapped Hilbert phase angles is usually
%   not ideal because phase is circular and has discontinuities at +/-pi.
%   To make the analysis more compatible with linear MVAR assumptions, this
%   script uses a linearized phase variable by default:
%
%       phaseVariable = sin(angle(hilbert(filteredSignal)))
%
%   and amplitude as:
%
%       amplitudeVariable = zscore(log(abs(hilbert(filteredSignal)) + eps))
%
%   This should be treated as an exploratory analysis, not a final
%   confirmatory one. A stronger confirmatory workflow would use
%   source-level signals and/or dedicated phase-amplitude coupling methods.
%
% Example:
%   results = computeCrossFrequencyPhaseAmplitudeGrangerBK1('013AR','M1','post');
%
%   cfg = struct;
%   cfg.modelOrder = 10;
%   cfg.saveOutput = 1;
%   results = computeCrossFrequencyPhaseAmplitudeGrangerBK1('013AR','M1','pre','ep','v8',cfg);
%
% Inputs:
%   subjectName      : BK1 subject code, e.g. '013AR'
%   protocolName     : one of EO1, EC1, G1, M1, G2, EO2, EC2, M2
%   conditionType    : 'pre' or 'post'
%   badEyeCondition  : default 'ep'
%   badTrialVersion  : default 'v8'
%   cfgIn            : optional struct
%       .ftDataFolder          default <root>/data/ftData
%       .stRange              default [0.25 1.25]
%       .bandEdges            default [4 8 12 16 32 36 40 44]
%       .targetBandPairs      default [1 5; 1 6; 1 7;
%                                       2 5; 2 6; 2 7;
%                                       3 5; 3 6; 3 7]
%                              Each row is [phaseBandIndex ampBandIndex].
%                              Indices refer to bands formed by bandEdges.
%                              Default bands: 1=4-8Hz, 2=8-12Hz, 3=12-16Hz,
%                              4=16-32Hz, 5=32-36Hz, 6=36-40Hz, 7=40-44Hz.
%       .filterOrder          default 4
%       .modelOrder           default 10
%       .saveOutput           default 1
%       .savedDataFolder      default <this folder>/savedDataCrossFreqGranger
%       .phaseRepresentation  default 'sine' ('sine' or 'raw')
%       .amplitudeTransform   default 'log'
%       .zscoreWithinTrial    default 1
%       .useParallel          default 1
%
% Output:
%   results structure with ROI labels, band definitions, and time-domain
%   forward Granger values for:
%       phase(source ROI, source band) -> amplitude(target ROI, target band)
%   Only the specified band pairs are computed (all 6x6 ROI combinations).

if ~exist('subjectName','var') || isempty(subjectName);         subjectName = '013AR'; end
if ~exist('protocolName','var') || isempty(protocolName);       protocolName = 'M1';    end
if ~exist('conditionType','var') || isempty(conditionType);     conditionType = 'post'; end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('cfgIn','var') || isempty(cfgIn);                     cfgIn = struct;          end

cfgIn = applyDefaultCfg(cfgIn);

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
ftDataFileName = fullfile(cfgIn.ftDataFolder,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '.mat']);
savedDataFolder = fullfile(cfgIn.savedDataFolder,subjectName);
ensureFolder(savedDataFolder);

if ~exist(ftDataFileName,'file')
    error('Could not find ftData file: %s', ftDataFileName);
end

tmp = load(ftDataFileName);
if ~isfield(tmp,'data') || isempty(tmp.data)
    error('No data field found in ftData file: %s', ftDataFileName);
end

data = tmp.data;
data = selectConditionWindow(data,conditionType,cfgIn.stRange);

[roiGroups,roiLabels] = getDefaultROIGroups();
bandInfo = makeBandInfoFromEdges(cfgIn.bandEdges);
targetBandPairs = cfgIn.targetBandPairs;

fs = data.fsample;
numROIs = numel(roiLabels);
numBands = numel(bandInfo);

roiChannelList = getValidROIChannels(data,roiGroups);
roiSignals = buildROISignals(data,roiChannelList);
phaseData = cell(numROIs,numBands);
ampData = cell(numROIs,numBands);
useParallel = shouldUseParallel(cfgIn);

uniqueAllBands = unique(targetBandPairs(:))';
if any(uniqueAllBands < 1) || any(uniqueAllBands > numBands)
    error('targetBandPairs contains band indices outside the valid range [1, %d].', numBands);
end

roiBandJobs = buildROIBandJobs(numROIs,uniqueAllBands);
phaseDataLinear = cell(size(roiBandJobs,1),1);
ampDataLinear = cell(size(roiBandJobs,1),1);

if useParallel
    parfor iJob = 1:size(roiBandJobs,1)
        iROI = roiBandJobs(iJob,1);
        iBand = roiBandJobs(iJob,2);
        [phaseDataLinear{iJob}, ampDataLinear{iJob}] = ...
            computeBandPhaseAmplitude(roiSignals{iROI}, bandInfo(iBand).range, fs, cfgIn);
    end
else
    for iJob = 1:size(roiBandJobs,1)
        iROI = roiBandJobs(iJob,1);
        iBand = roiBandJobs(iJob,2);
        [phaseDataLinear{iJob}, ampDataLinear{iJob}] = ...
            computeBandPhaseAmplitude(roiSignals{iROI}, bandInfo(iBand).range, fs, cfgIn);
    end
end

for iJob = 1:size(roiBandJobs,1)
    iROI = roiBandJobs(iJob,1);
    iBand = roiBandJobs(iJob,2);
    phaseData{iROI,iBand} = phaseDataLinear{iJob};
    ampData{iROI,iBand} = ampDataLinear{iJob};
end

gcJobs = buildROIBandPairJobs(numROIs,targetBandPairs);
numGCJobs = size(gcJobs,1);
phaseToAmpGCLinear = nan(numGCJobs,1);
sampleCountLinear = nan(numGCJobs,1);
residualVarFullLinear = nan(numGCJobs,2);
residualVarRestrictedLinear = nan(numGCJobs,1);

if useParallel
    parfor iJob = 1:numGCJobs
        iSrcROI = gcJobs(iJob,1);
        iDstROI = gcJobs(iJob,2);
        iPhaseBand = gcJobs(iJob,3);
        iAmpBand = gcJobs(iJob,4);

        phaseTrials = phaseData{iSrcROI,iPhaseBand};
        ampTrials = ampData{iDstROI,iAmpBand};

        if ~isempty(phaseTrials) && ~isempty(ampTrials)
            [gcForward,nObs,fullVar,restVarY] = ...
                computeForwardTimeDomainGC(phaseTrials,ampTrials,cfgIn.modelOrder);
            phaseToAmpGCLinear(iJob) = gcForward;
            sampleCountLinear(iJob) = nObs;
            residualVarFullLinear(iJob,:) = fullVar;
            residualVarRestrictedLinear(iJob) = restVarY;
        end
    end
else
    for iJob = 1:numGCJobs
        iSrcROI = gcJobs(iJob,1);
        iDstROI = gcJobs(iJob,2);
        iPhaseBand = gcJobs(iJob,3);
        iAmpBand = gcJobs(iJob,4);

        phaseTrials = phaseData{iSrcROI,iPhaseBand};
        ampTrials = ampData{iDstROI,iAmpBand};

        if ~isempty(phaseTrials) && ~isempty(ampTrials)
            [gcForward,nObs,fullVar,restVarY] = ...
                computeForwardTimeDomainGC(phaseTrials,ampTrials,cfgIn.modelOrder);
            phaseToAmpGCLinear(iJob) = gcForward;
            sampleCountLinear(iJob) = nObs;
            residualVarFullLinear(iJob,:) = fullVar;
            residualVarRestrictedLinear(iJob) = restVarY;
        end
    end
end

phaseToAmpGC = nan(numROIs,numROIs,numBands,numBands);
sampleCount = nan(numROIs,numROIs,numBands,numBands);
residualVarFull = nan(numROIs,numROIs,numBands,numBands,2);
residualVarRestricted = nan(numROIs,numROIs,numBands,numBands,2);

for iJob = 1:numGCJobs
    iSrcROI = gcJobs(iJob,1);
    iDstROI = gcJobs(iJob,2);
    iPhaseBand = gcJobs(iJob,3);
    iAmpBand = gcJobs(iJob,4);

    phaseToAmpGC(iSrcROI,iDstROI,iPhaseBand,iAmpBand) = phaseToAmpGCLinear(iJob);
    sampleCount(iSrcROI,iDstROI,iPhaseBand,iAmpBand) = sampleCountLinear(iJob);
    residualVarFull(iSrcROI,iDstROI,iPhaseBand,iAmpBand,:) = residualVarFullLinear(iJob,:);
    residualVarRestricted(iSrcROI,iDstROI,iPhaseBand,iAmpBand,1) = residualVarRestrictedLinear(iJob);
end

results = struct;
results.subjectName = subjectName;
results.protocolName = protocolName;
results.conditionType = conditionType;
results.badEyeCondition = badEyeCondition;
results.badTrialVersion = badTrialVersion;
results.numGoodTrials = tmp.numGoodTrials;
results.fs = fs;
results.stRange = cfgIn.stRange;
results.phaseRepresentation = cfgIn.phaseRepresentation;
results.amplitudeTransform = cfgIn.amplitudeTransform;
results.modelOrder = cfgIn.modelOrder;
results.bandInfo = bandInfo;
results.targetBandPairs = targetBandPairs;
results.roiLabels = roiLabels;
results.roiGroups = roiGroups;
results.roiChannelList = roiChannelList;
results.phaseToAmpGC = phaseToAmpGC;
results.sampleCount = sampleCount;
results.residualVarFull = residualVarFull;
results.residualVarRestricted = residualVarRestricted;
results.notes = ['Exploratory regional phase-amplitude Granger. ', ...
                 'Forward direction is phase(source ROI/source band) -> amplitude(target ROI/target band). ', ...
                 'Only the 9 pre-specified low-frequency phase / high-frequency amplitude band pairs are computed. ', ...
                 'Phase variable is linearized for MVAR compatibility. ', ...
                 'residualVarRestricted slot 1 = restricted Y variance. ', ...
                 'residualVarRestricted slot 2 = NaN because reverse GC is not computed.'];

if cfgIn.saveOutput
    saveFileName = fullfile(savedDataFolder,[protocolName '_' conditionType '_' badEyeCondition '_' badTrialVersion '_crossFreqPhaseAmpGC.mat']);
    save(saveFileName,'results','-v7.3');
    fprintf('Saved cross-frequency phase-amplitude GC to:\n%s\n',saveFileName);
end
end

function cfgOut = applyDefaultCfg(cfgIn)
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));

cfgOut = cfgIn;
cfgOut = setDefault(cfgOut,'ftDataFolder',fullfile(projectRoot,'data','ftData'));
cfgOut = setDefault(cfgOut,'stRange',[0.25 1.25]);
cfgOut = setDefault(cfgOut,'bandEdges',[4 8 12 16 32 36 40 44]);
cfgOut = setDefault(cfgOut,'targetBandPairs', ...
    [1 5; ...
     1 6; ...
     1 7; ...
     2 5; ...
     2 6; ...
     2 7; ...
     3 5; ...
     3 6; ...
     3 7]);
cfgOut = setDefault(cfgOut,'filterOrder',4);
cfgOut = setDefault(cfgOut,'modelOrder',10);
cfgOut = setDefault(cfgOut,'saveOutput',1);
cfgOut = setDefault(cfgOut,'savedDataFolder',fullfile(fileparts(mfilename('fullpath')),'savedDataCrossFreqGranger'));
cfgOut = setDefault(cfgOut,'phaseRepresentation','sine');
cfgOut = setDefault(cfgOut,'amplitudeTransform','log');
cfgOut = setDefault(cfgOut,'zscoreWithinTrial',1);
cfgOut = setDefault(cfgOut,'useParallel',1);
end

function s = setDefault(s,fieldName,defaultValue)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end

function dataCond = selectConditionWindow(data,conditionType,stRange)
cfg = [];
switch lower(conditionType)
    case 'pre'
        cfg.toilim = [(-diff(stRange) + 1/data.fsample) 0];
    case 'post'
        cfg.toilim = [(stRange(1) + 1/data.fsample) stRange(2)];
    otherwise
        error('conditionType must be ''pre'' or ''post''.');
end
dataCond = ft_redefinetrial(cfg,data);
end

function [roiGroups,roiLabels] = getDefaultROIGroups()
roiGroups{1} = [14:16 32+(12:15)];             roiLabels{1} = 'Occipital_L';
roiGroups{2} = [18:20 32+(17:20)];             roiLabels{2} = 'Occipital_R';
roiGroups{3} = [6:8 11 12 32+[7:9 11]];        roiLabels{3} = 'Central_L';
roiGroups{4} = [22 23 25 28 29 32+[22 24:26]]; roiLabels{4} = 'Central_R';
roiGroups{5} = [1 3 4 32+[1 2 4 5]];           roiLabels{5} = 'Frontal_L';
roiGroups{6} = [30:32 32+(28:31)];             roiLabels{6} = 'Frontal_R';
end

function bandInfo = makeBandInfoFromEdges(bandEdges)
bandEdges = bandEdges(:)';
if numel(bandEdges) < 2
    error('bandEdges must contain at least two values.');
end

numBands = numel(bandEdges) - 1;
bandInfo = struct('name',cell(1,numBands),'range',cell(1,numBands));
for iBand = 1:numBands
    bandInfo(iBand).range = [bandEdges(iBand) bandEdges(iBand+1)];
    bandInfo(iBand).name = sprintf('%g_%gHz',bandEdges(iBand),bandEdges(iBand+1));
end
end

function roiChannelList = getValidROIChannels(data,roiGroups)
numAllChannels = numel(data.label);
badElecs = [];
if isfield(data,'badElecs') && ~isempty(data.badElecs)
    badElecs = unique(data.badElecs(:)');
    badElecs = badElecs(badElecs >= 1 & badElecs <= numAllChannels);
end

roiChannelList = cell(size(roiGroups));
for iROI = 1:numel(roiGroups)
    roiChannelList{iROI} = setdiff(roiGroups{iROI},badElecs,'stable');
end
end

function roiSignals = buildROISignals(data,roiChannelList)
numROIs = numel(roiChannelList);
numTrials = numel(data.trial);
roiSignals = cell(1,numROIs);

for iROI = 1:numROIs
    roiSignals{iROI} = cell(1,numTrials);
    for iTrial = 1:numTrials
        trialData = data.trial{iTrial};
        useChans = roiChannelList{iROI};
        if isempty(useChans)
            roiSignals{iROI}{iTrial} = [];
        else
            roiSignals{iROI}{iTrial} = mean(trialData(useChans,:),1,'omitnan');
        end
    end
end
end

function [phaseTrials, ampTrials] = computeBandPhaseAmplitude(signalTrials,bandRange,fs,cfgIn)
numTrials = numel(signalTrials);
phaseTrials = cell(1,numTrials);
ampTrials = cell(1,numTrials);

for iTrial = 1:numTrials
    x = signalTrials{iTrial};
    if isempty(x) || all(~isfinite(x))
        phaseTrials{iTrial} = [];
        ampTrials{iTrial} = [];
        continue;
    end

    x = x(:)';
    x = x - mean(x,'omitnan');
    xf = bandpassFilterTrial(x,bandRange,fs,cfgIn.filterOrder);
    analyticSig = hilbert(xf);

    rawPhase = angle(analyticSig);
    switch lower(cfgIn.phaseRepresentation)
        case 'sine'
            phaseVals = sin(rawPhase);
        case 'raw'
            phaseVals = rawPhase;
        otherwise
            error('Unknown phaseRepresentation: %s',cfgIn.phaseRepresentation);
    end

    ampVals = abs(analyticSig);
    switch lower(cfgIn.amplitudeTransform)
        case 'log'
            ampVals = log(ampVals + eps);
        case 'none'
            % do nothing
        otherwise
            error('Unknown amplitudeTransform: %s',cfgIn.amplitudeTransform);
    end

    if cfgIn.zscoreWithinTrial
        phaseVals = safeZscore(phaseVals);
        ampVals = safeZscore(ampVals);
    end

    phaseTrials{iTrial} = phaseVals(:)';
    ampTrials{iTrial} = ampVals(:)';
end
end

function jobList = buildROIBandJobs(numROIs,bandIndices)
[roiIdx,bandIdx] = ndgrid(1:numROIs,bandIndices);
jobList = [roiIdx(:) bandIdx(:)];
end

function jobList = buildROIBandPairJobs(numROIs,targetBandPairs)
numPairs = size(targetBandPairs,1);
jobList = nan(numROIs*numROIs*numPairs,4);
jobCounter = 0;
for iSrcROI = 1:numROIs
    for iDstROI = 1:numROIs
        for iPair = 1:numPairs
            jobCounter = jobCounter + 1;
            jobList(jobCounter,:) = [iSrcROI iDstROI targetBandPairs(iPair,1) targetBandPairs(iPair,2)];
        end
    end
end
end

function tf = shouldUseParallel(cfgIn)
tf = false;
if ~isfield(cfgIn,'useParallel') || ~cfgIn.useParallel
    return;
end
if exist('gcp','file') ~= 2 || exist('parfor','builtin') ~= 5
    return;
end
try
    tf = ~isempty(gcp('nocreate'));
catch
    tf = false;
end
end

function xf = bandpassFilterTrial(x,bandRange,fs,filterOrder)
wn = bandRange ./ (fs/2);
wn(1) = max(wn(1),1e-6);
wn(2) = min(wn(2),0.999);
if wn(1) >= wn(2)
    error('Invalid normalized band range.');
end
[b,a] = butter(filterOrder,wn,'bandpass');
xf = filtfilt(b,a,double(x));
end

function xz = safeZscore(x)
mu = mean(x,'omitnan');
sd = std(x,0,'omitnan');
if ~isfinite(sd) || sd == 0
    xz = x - mu;
else
    xz = (x - mu) ./ sd;
end
end

function [gcXtoY,nObs,fullVar,varRestY] = computeForwardTimeDomainGC(xTrials,yTrials,modelOrder)
[Xfull,Yfull] = buildFullModelRows(xTrials,yTrials,modelOrder);
[XrestY,Yresp] = buildRestrictedModelRows(xTrials,yTrials,modelOrder,'y');

nObs = size(Yfull,1);
gcXtoY = NaN;
fullVar = [NaN NaN];
varRestY = NaN;

if nObs <= (2*modelOrder + 1)
    return;
end

betaFull = Xfull \ Yfull;
resFull = Yfull - Xfull*betaFull;
varFullX = var(resFull(:,1),1);
varFullY = var(resFull(:,2),1);

betaRestY = XrestY \ Yresp;
resRestY = Yresp - XrestY*betaRestY;
varRestY = var(resRestY,1);

if isfinite(varFullY) && isfinite(varRestY) && varFullY > 0 && varRestY > 0
    gcXtoY = log(varRestY / varFullY);
end

fullVar = [varFullX varFullY];
end

function [Xfull,Yfull] = buildFullModelRows(xTrials,yTrials,p)
Xfull = [];
Yfull = [];

for iTrial = 1:numel(xTrials)
    x = xTrials{iTrial};
    y = yTrials{iTrial};
    if isempty(x) || isempty(y)
        continue;
    end
    T = min(numel(x),numel(y));
    x = x(1:T);
    y = y(1:T);
    if T <= p
        continue;
    end

    Xt = nan(T-p,2*p + 1);
    Yt = nan(T-p,2);

    for t = (p+1):T
        row = t - p;
        xLags = x(t-1:-1:t-p);
        yLags = y(t-1:-1:t-p);
        Xt(row,:) = [1 xLags yLags];
        Yt(row,:) = [x(t) y(t)];
    end

    validRows = all(isfinite(Xt),2) & all(isfinite(Yt),2);
    Xfull = [Xfull; Xt(validRows,:)]; %#ok<AGROW>
    Yfull = [Yfull; Yt(validRows,:)]; %#ok<AGROW>
end
end

function [Xrest,response] = buildRestrictedModelRows(xTrials,yTrials,p,responseType)
Xrest = [];
response = [];

for iTrial = 1:numel(xTrials)
    x = xTrials{iTrial};
    y = yTrials{iTrial};
    if isempty(x) || isempty(y)
        continue;
    end
    T = min(numel(x),numel(y));
    x = x(1:T);
    y = y(1:T);
    if T <= p
        continue;
    end

    Xt = nan(T-p,p + 1);
    Rt = nan(T-p,1);

    for t = (p+1):T
        row = t - p;
        xLags = x(t-1:-1:t-p);
        yLags = y(t-1:-1:t-p);

        switch lower(responseType)
            case 'y'
                Xt(row,:) = [1 yLags];
                Rt(row) = y(t);
            case 'x'
                Xt(row,:) = [1 xLags];
                Rt(row) = x(t);
            otherwise
                error('responseType must be ''x'' or ''y''.');
        end
    end

    validRows = all(isfinite(Xt),2) & isfinite(Rt);
    Xrest = [Xrest; Xt(validRows,:)]; %#ok<AGROW>
    response = [response; Rt(validRows)]; %#ok<AGROW>
end
end

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
