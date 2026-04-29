function results = extractGrangerBandDataAllSubjectsBK1(protocolName,badEyeCondition,badTrialVersion,saveOutputFlag,bandInfo)
% extractGrangerBandDataAllSubjectsBK1 Extract band-averaged BK1 Granger data.
%
% Example:
%   results = extractGrangerBandDataAllSubjectsBK1('M1');
%   bandInfo(1).name = 'custom1'; bandInfo(1).range = [25 45];
%   results = extractGrangerBandDataAllSubjectsBK1('M1','ep','v8',1,bandInfo);
%
% Output:
%   Saves a summary MAT file in savedDataGranger/summary containing
%   channel-level and ROI-level pre, post, and post-pre Granger matrices.

if ~exist('protocolName','var') || isempty(protocolName);       protocolName = 'M1';    end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('saveOutputFlag','var') || isempty(saveOutputFlag);   saveOutputFlag = 1;     end
if ~exist('bandInfo','var') || isempty(bandInfo);               bandInfo = getDefaultBands(); end

projectRoot = getProjectRoot();
grangerFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataGranger');
summaryFolder = fullfile(grangerFolder,'summary');
ensureFolder(summaryFolder);

[subjectNameList,groupNameList] = getBK1SubjectGroups(projectRoot);
[roiGroups,roiLabels] = getDefaultROIGroups();
bandInfo = validateBandInfo(bandInfo);
electrodeLabels = getDefaultElectrodeLabels(projectRoot);
bandTag = getBandInfoTag(bandInfo);

numSubjects = numel(subjectNameList);
numElectrodes = numel(electrodeLabels);
numROIs = numel(roiLabels);
numBands = numel(bandInfo);

results = struct;
results.protocolName = protocolName;
results.badEyeCondition = badEyeCondition;
results.badTrialVersion = badTrialVersion;
results.sourceFolder = grangerFolder;
results.summaryFolder = summaryFolder;
results.subjectNameList = subjectNameList;
results.groupNameList = groupNameList;
results.electrodeLabels = electrodeLabels;
results.roiLabels = roiLabels;
results.roiGroups = roiGroups;
results.bandInfo = bandInfo;
results.numGoodTrials = nan(numSubjects,1);
results.samplesPerParameterPre = nan(numSubjects,1);
results.samplesPerParameterPost = nan(numSubjects,1);
results.lowSampleRatioFlagPre = nan(numSubjects,1);
results.lowSampleRatioFlagPost = nan(numSubjects,1);
results.validSubjectMask = false(numSubjects,1);
results.preStatusList = cell(numSubjects,1);
results.postStatusList = cell(numSubjects,1);
results.selectedElectrodeList = cell(numSubjects,1);
results.bandTag = bandTag;
results.outputFileName = fullfile(summaryFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_' bandTag '_grangerBandData.mat']);

for iBand = 1:numBands
    results.bandData(iBand).name = bandInfo(iBand).name;
    results.bandData(iBand).range = bandInfo(iBand).range;
    results.bandData(iBand).channelPre = nan(numElectrodes,numElectrodes,numSubjects);
    results.bandData(iBand).channelPost = nan(numElectrodes,numElectrodes,numSubjects);
    results.bandData(iBand).channelDiff = nan(numElectrodes,numElectrodes,numSubjects);
    results.bandData(iBand).roiPre = nan(numROIs,numROIs,numSubjects);
    results.bandData(iBand).roiPost = nan(numROIs,numROIs,numSubjects);
    results.bandData(iBand).roiDiff = nan(numROIs,numROIs,numSubjects);
    results.bandData(iBand).channelNetPre = nan(numElectrodes,numSubjects);
    results.bandData(iBand).channelNetPost = nan(numElectrodes,numSubjects);
    results.bandData(iBand).channelNetDiff = nan(numElectrodes,numSubjects);
    results.bandData(iBand).roiNetPre = nan(numROIs,numSubjects);
    results.bandData(iBand).roiNetPost = nan(numROIs,numSubjects);
    results.bandData(iBand).roiNetDiff = nan(numROIs,numSubjects);
end

for iSubject = 1:numSubjects
    subjectName = subjectNameList{iSubject};
    dataFileName = fullfile(grangerFolder,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '_granger.mat']);

    if ~exist(dataFileName,'file')
        results.preStatusList{iSubject} = 'missing_file';
        results.postStatusList{iSubject} = 'missing_file';
        continue;
    end

    x = load(dataFileName);
    results.numGoodTrials(iSubject) = x.numGoodTrials;
    results.preStatusList{iSubject} = getStatusString(x.metaPre);
    results.postStatusList{iSubject} = getStatusString(x.metaPost);
    results.samplesPerParameterPre(iSubject) = getMetaField(x.metaPre,'samplesPerParameter');
    results.samplesPerParameterPost(iSubject) = getMetaField(x.metaPost,'samplesPerParameter');
    results.lowSampleRatioFlagPre(iSubject) = getMetaField(x.metaPre,'lowSampleRatioFlag');
    results.lowSampleRatioFlagPost(iSubject) = getMetaField(x.metaPost,'lowSampleRatioFlag');
    results.selectedElectrodeList{iSubject} = intersect(getSelectedElectrodes(x.metaPre),getSelectedElectrodes(x.metaPost),'stable');

    if ~strcmp(results.preStatusList{iSubject},'ok') || ~strcmp(results.postStatusList{iSubject},'ok')
        continue;
    end

    results.validSubjectMask(iSubject) = true;

    for iBand = 1:numBands
        preMask = x.freqPre >= bandInfo(iBand).range(1) & x.freqPre <= bandInfo(iBand).range(2);
        postMask = x.freqPost >= bandInfo(iBand).range(1) & x.freqPost <= bandInfo(iBand).range(2);

        if ~any(preMask) || ~any(postMask)
            continue;
        end

        Gpre = averageBandMatrix(x.grangerPre,preMask);
        Gpost = averageBandMatrix(x.grangerPost,postMask);
        Gdiff = Gpost - Gpre;

        Gpre = removeSelfConnections(Gpre);
        Gpost = removeSelfConnections(Gpost);
        Gdiff = removeSelfConnections(Gdiff);

        results.bandData(iBand).channelPre(:,:,iSubject) = Gpre;
        results.bandData(iBand).channelPost(:,:,iSubject) = Gpost;
        results.bandData(iBand).channelDiff(:,:,iSubject) = Gdiff;

        results.bandData(iBand).roiPre(:,:,iSubject) = computeROIMatrix(Gpre,roiGroups);
        results.bandData(iBand).roiPost(:,:,iSubject) = computeROIMatrix(Gpost,roiGroups);
        results.bandData(iBand).roiDiff(:,:,iSubject) = computeROIMatrix(Gdiff,roiGroups);

        results.bandData(iBand).channelNetPre(:,iSubject) = computeNetFlow(Gpre);
        results.bandData(iBand).channelNetPost(:,iSubject) = computeNetFlow(Gpost);
        results.bandData(iBand).channelNetDiff(:,iSubject) = computeNetFlow(Gdiff);

        results.bandData(iBand).roiNetPre(:,iSubject) = computeNetFlow(results.bandData(iBand).roiPre(:,:,iSubject));
        results.bandData(iBand).roiNetPost(:,iSubject) = computeNetFlow(results.bandData(iBand).roiPost(:,:,iSubject));
        results.bandData(iBand).roiNetDiff(:,iSubject) = computeNetFlow(results.bandData(iBand).roiDiff(:,:,iSubject));
    end
end

if saveOutputFlag
    save(results.outputFileName,'results','-v7.3');
end

disp(sprintf('Extracted Granger band data for protocol %s: %d/%d valid subjects.',protocolName,sum(results.validSubjectMask),numSubjects)); %#ok<DSPS>
end

function projectRoot = getProjectRoot()
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

function [subjectNameList,groupNameList] = getBK1SubjectGroups(projectRoot)
infoFileName = fullfile(projectRoot,'ProjectDhyaanBK1Programs','commonAnalysisCodes','informationFiles','BK1AllSubjectList.mat');
tmp = load(infoFileName);

declaredBadSubjects = {'004P','081SN','069MG','092KB','063VK'};
allSubjectList = setdiff(normalizeList(tmp.allSubjectList),declaredBadSubjects,'stable');
meditatorList = setdiff(normalizeList(tmp.meditatorList),declaredBadSubjects,'stable');
controlList = setdiff(normalizeList(tmp.controlList),declaredBadSubjects,'stable');

subjectNameList = allSubjectList(:);
groupNameList = cell(size(subjectNameList));
for iSubject = 1:numel(subjectNameList)
    if ismember(subjectNameList{iSubject},meditatorList)
        groupNameList{iSubject} = 'Meditator';
    elseif ismember(subjectNameList{iSubject},controlList)
        groupNameList{iSubject} = 'Control';
    else
        groupNameList{iSubject} = 'Unknown';
    end
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

function electrodeLabels = getDefaultElectrodeLabels(projectRoot)
labelFileName = fullfile(projectRoot,'Montages','Layouts','actiCap64_UOL','actiCap64_UOLLabels.mat');
tmp = load(labelFileName);
electrodeLabels = tmp.montageLabels(:,2);
end

function [roiGroups,roiLabels] = getDefaultROIGroups()
roiGroups{1} = [14:16 32+(12:15)];         roiLabels{1} = 'Occipital_L';
roiGroups{2} = [18:20 32+(17:20)];         roiLabels{2} = 'Occipital_R';
roiGroups{3} = [6:8 11 12 32+[7:9 11]];    roiLabels{3} = 'Central_L';
roiGroups{4} = [22 23 25 28 29 32+[22 24:26]]; roiLabels{4} = 'Central_R';
roiGroups{5} = [1 3 4 32+[1 2 4 5]];       roiLabels{5} = 'Frontal_L';
roiGroups{6} = [30:32 32+(28:31)];         roiLabels{6} = 'Frontal_R';
end

function bandInfo = getDefaultBands()
bandInfo(1).name = 'alpha';     bandInfo(1).range = [7 10];
bandInfo(2).name = 'beta';      bandInfo(2).range = [20 32];
bandInfo(3).name = 'highgamma'; bandInfo(3).range = [30 80];
end

function bandInfo = validateBandInfo(bandInfo)
if ~isstruct(bandInfo) || isempty(bandInfo)
    error('bandInfo must be a non-empty struct array with fields name and range.');
end

for iBand = 1:numel(bandInfo)
    if ~isfield(bandInfo(iBand),'name') || ~isfield(bandInfo(iBand),'range')
        error('Each bandInfo entry must contain name and range.');
    end
    if numel(bandInfo(iBand).range) ~= 2
        error('Each bandInfo range must have exactly two values.');
    end
    bandInfo(iBand).range = reshape(bandInfo(iBand).range,1,2);
    if bandInfo(iBand).range(1) >= bandInfo(iBand).range(2)
        error('Band range must satisfy low < high.');
    end
    bandInfo(iBand).name = char(string(bandInfo(iBand).name));
end
end

function bandTag = getBandInfoTag(bandInfo)
tagParts = cell(1,numel(bandInfo));
for iBand = 1:numel(bandInfo)
    safeName = regexprep(lower(bandInfo(iBand).name),'[^a-z0-9]+','');
    rangeVals = round(100*bandInfo(iBand).range);
    tagParts{iBand} = sprintf('%s_%g_%g',safeName,rangeVals(1),rangeVals(2));
end
bandTag = strjoin(tagParts,'__');
end

function statusString = getStatusString(metaStruct)
statusString = 'missing_meta';
if isfield(metaStruct,'status')
    statusString = metaStruct.status;
end
end

function val = getMetaField(metaStruct,fieldName)
val = NaN;
if isfield(metaStruct,fieldName)
    val = metaStruct.(fieldName);
end
end

function selectedElectrodes = getSelectedElectrodes(metaStruct)
selectedElectrodes = [];
if isfield(metaStruct,'selectedElectrodes')
    selectedElectrodes = metaStruct.selectedElectrodes(:)';
end
end

function G = averageBandMatrix(grangerData,freqMask)
G = mean(grangerData(:,:,freqMask),3,'omitnan');
end

function G = removeSelfConnections(G)
n = size(G,1);
G(1:n+1:end) = NaN;
end

function roiMatrix = computeROIMatrix(channelMatrix,roiGroups)
numROIs = numel(roiGroups);
roiMatrix = nan(numROIs,numROIs);

for iROI = 1:numROIs
    srcIdx = roiGroups{iROI};
    for jROI = 1:numROIs
        dstIdx = roiGroups{jROI};
        block = channelMatrix(srcIdx,dstIdx);
        if iROI == jROI
            block(1:size(block,1)+1:end) = NaN;
        end
        roiMatrix(iROI,jROI) = mean(block(:),'omitnan');
    end
end
end

function netFlow = computeNetFlow(channelMatrix)
outFlow = mean(channelMatrix,2,'omitnan');
inFlow = mean(channelMatrix,1,'omitnan')';
netFlow = outFlow - inFlow;
end

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
