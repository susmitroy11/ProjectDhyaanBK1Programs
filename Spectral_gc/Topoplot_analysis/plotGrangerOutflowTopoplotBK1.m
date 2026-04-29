function topo = plotGrangerOutflowTopoplotBK1(selection,protocolName,freqRange,conditionType,badEyeCondition,badTrialVersion,saveFigureFlag,useMedianFlag,measureType)
% plotGrangerOutflowTopoplotBK1 Plot BK1 scalp topographies of Granger flow.
% This function can summarize outflow, inflow, or net flow across a chosen
% frequency range for one subject, one group, or a batch of protocols and
% conditions.
%
% Examples:
%   plotGrangerOutflowTopoplotBK1('028HB','M1',[30 80],'diff');
%   plotGrangerOutflowTopoplotBK1('Meditator','M1',[30 80],'diff');
%   plotGrangerOutflowTopoplotBK1('Control','M1',[30 80],'post', 'ep', 'v8', 1);
%   plotGrangerOutflowTopoplotBK1('Meditator','M1',[30 80],'diff', 'ep', 'v8', 1, 0, 'netflow');
%   plotGrangerOutflowTopoplotBK1('BothGroups','all',[30 80],'all');
%
% selection:
%   subject code like '028HB'
%   or one of: 'Meditator', 'Control', 'All', 'BothGroups'
%
% conditionType:
%   'pre'  = average pre-window outflow
%   'post' = average post-window outflow
%   'diff' = post - pre outflow (default)
%   'all'  = generate both pre and post plots
%
% protocolName:
%   protocol code like 'M1'
%   or 'all' to generate all BK1 protocols

if ~exist('selection','var') || isempty(selection);            selection = 'Meditator'; end
if ~exist('protocolName','var') || isempty(protocolName);      protocolName = 'M1';     end
if ~exist('freqRange','var') || isempty(freqRange);            freqRange = [30 80];     end
if ~exist('conditionType','var') || isempty(conditionType);    conditionType = 'diff';  end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('saveFigureFlag','var') || isempty(saveFigureFlag);  saveFigureFlag = 0;      end
if ~exist('useMedianFlag','var') || isempty(useMedianFlag);    useMedianFlag = 0;       end
if ~exist('measureType','var') || isempty(measureType);        measureType = 'outflow'; end

if shouldRunBatchMode(selection,protocolName,conditionType)
    topo = runBatchTopoplot(selection,protocolName,freqRange,conditionType,badEyeCondition,badTrialVersion,saveFigureFlag,useMedianFlag,measureType);
    return;
end

projectRoot = getProjectRoot();
grangerFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataGranger');
figureFolder = fullfile(grangerFolder,'figures','flowTopoplots');
ensureFolder(figureFolder);

addEEGLABTopoplot(projectRoot);
[montageChanlocs,montageLabels] = getMontageChanlocs(projectRoot);
[subjectNameList,groupNameList] = getBK1SubjectGroups(projectRoot);
selectedSubjects = resolveSelection(selection,subjectNameList,groupNameList);

if isempty(selectedSubjects)
    error('No subjects found for selection: %s',char(selection));
end

numSubjects = numel(selectedSubjects);
numElectrodes = numel(montageLabels);
outflowAll = nan(numElectrodes,numSubjects);
inflowAll = nan(numElectrodes,numSubjects);
statusList = cell(numSubjects,1);

for iSubject = 1:numSubjects
    subjectName = selectedSubjects{iSubject};
    dataFileName = fullfile(grangerFolder,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '_granger.mat']);

    if ~exist(dataFileName,'file')
        statusList{iSubject} = 'missing_file';
        continue;
    end

    x = load(dataFileName);
    [outflow,inflow,statusString] = computeOutflowFromFile(x,freqRange,conditionType,numElectrodes);
    outflowAll(:,iSubject) = outflow;
    inflowAll(:,iSubject) = inflow;
    statusList{iSubject} = statusString;
end

validMask = strcmp(statusList,'ok');
if ~any(validMask)
    error('No valid Granger files found for the requested topoplot.');
end

outflowAll = outflowAll(:,validMask);
inflowAll = inflowAll(:,validMask);
usedSubjects = selectedSubjects(validMask);

if useMedianFlag
    topoValues = summarizeMeasure(outflowAll,inflowAll,measureType,1);
    summaryMethod = 'median';
else
    topoValues = summarizeMeasure(outflowAll,inflowAll,measureType,0);
    summaryMethod = 'mean';
end
sharedColorLimits = [-0.07 0.07];

[maxValue,maxPositiveIndex] = max(topoValues);
[minValue,maxNegativeIndex] = min(topoValues);

hFig = figure('Color','w','Name',['Granger ' capitalizeLabel(measureType) ' Topoplot']);
hFig.Position(3:4) = [900 700];
topoplot(topoValues,montageChanlocs,'electrodes','on','plotrad',0.6,'headrad',0.6);
caxis(sharedColorLimits);
colorbar;
title(makeTitle(selection,protocolName,freqRange,conditionType,summaryMethod,numel(usedSubjects),measureType),'Interpreter','none');

topo = struct;
topo.selection = selection;
topo.protocolName = protocolName;
topo.freqRange = freqRange;
topo.conditionType = conditionType;
topo.summaryMethod = summaryMethod;
topo.subjectNameList = usedSubjects;
topo.statusList = statusList(validMask);
topo.outflowAll = outflowAll;
topo.inflowAll = inflowAll;
topo.measureType = measureType;
topo.topoValues = topoValues;
topo.colorLimits = sharedColorLimits;
topo.maxPositiveElectrode = montageLabels{maxPositiveIndex};
topo.maxPositiveIndex = maxPositiveIndex;
topo.maxPositiveValue = maxValue;
topo.maxNegativeElectrode = montageLabels{maxNegativeIndex};
topo.maxNegativeIndex = maxNegativeIndex;
topo.maxNegativeValue = minValue;
topo.figureHandle = hFig;

disp(sprintf('Top positive %s electrode: %s',measureType,topo.maxPositiveElectrode)); %#ok<DSPS>
disp(sprintf('Top negative %s electrode: %s',measureType,topo.maxNegativeElectrode)); %#ok<DSPS>

if saveFigureFlag
    protocolFigureFolder = fullfile(figureFolder,protocolName);
    ensureFolder(protocolFigureFolder);
    fileTag = sanitizeFileTag(selection);
    saveFileName = fullfile(protocolFigureFolder,[protocolName '_' lower(measureType) '_' conditionType '_' fileTag '_' num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz_topoplot.png']);
    saveas(hFig,saveFileName);
    topo.saveFileName = saveFileName;
end
end

function topoList = runBatchTopoplot(selection,protocolName,freqRange,conditionType,badEyeCondition,badTrialVersion,saveFigureFlag,useMedianFlag,measureType)
selectionList = expandBatchSelections(selection);
protocolList = expandBatchProtocols(protocolName);
conditionList = expandBatchConditions(conditionType);

numPlots = numel(selectionList) * numel(protocolList) * numel(conditionList);
topoList = cell(1,numPlots);
plotCounter = 0;

for iProtocol = 1:numel(protocolList)
    for iCondition = 1:numel(conditionList)
        for iSelection = 1:numel(selectionList)
            plotCounter = plotCounter + 1;
            topoList{plotCounter} = plotGrangerOutflowTopoplotBK1(selectionList{iSelection},protocolList{iProtocol},freqRange,conditionList{iCondition},badEyeCondition,badTrialVersion,saveFigureFlag,useMedianFlag,measureType);
        end
    end
end

disp(sprintf('Generated %d %s topoplots across %d protocols, %d conditions, and %d group selections.', ... %#ok<DSPS>
    numel(topoList),measureType,numel(protocolList),numel(conditionList),numel(selectionList)));
end

function tf = shouldRunBatchMode(selection,protocolName,conditionType)
tf = isBatchSelection(selection) || strcmpi(char(protocolName),'all') || strcmpi(char(conditionType),'all');
end

function tf = isBatchSelection(selection)
selectionText = char(selection);
tf = any(strcmpi(selectionText,{'BothGroups','Groups','Grouped'}));
end

function selectionList = expandBatchSelections(selection)
if isBatchSelection(selection)
    selectionList = {'Meditator','Control'};
else
    selectionList = {char(selection)};
end
end

function protocolList = expandBatchProtocols(protocolName)
if strcmpi(char(protocolName),'all')
    protocolList = getDefaultProtocolList();
else
    protocolList = {char(protocolName)};
end
end

function conditionList = expandBatchConditions(conditionType)
if strcmpi(char(conditionType),'all')
    conditionList = {'pre','post'};
else
    conditionList = {char(conditionType)};
end
end

function [outflow,inflow,statusString] = computeOutflowFromFile(x,freqRange,conditionType,numElectrodes)
statusString = 'ok';
outflow = nan(numElectrodes,1);
inflow = nan(numElectrodes,1);

preMask = x.freqPre >= freqRange(1) & x.freqPre <= freqRange(2);
postMask = x.freqPost >= freqRange(1) & x.freqPost <= freqRange(2);
preStatusOk = strcmp(getStatusString(x.metaPre),'ok');
postStatusOk = strcmp(getStatusString(x.metaPost),'ok');

switch lower(conditionType)
    case 'pre'
        if ~preStatusOk
            statusString = 'invalid_pre_meta';
            return;
        end
        if ~any(preMask)
            statusString = 'pre_freq_not_found';
            return;
        end
        Gpre = mean(x.grangerPre(:,:,preMask),3,'omitnan');
        Gpre = removeSelfConnections(Gpre);
        outflowPre = mean(Gpre,2,'omitnan');
        inflowPre = mean(Gpre,1,'omitnan')';
        outflow = outflowPre;
        inflow = inflowPre;
    case 'post'
        if ~postStatusOk
            statusString = 'invalid_post_meta';
            return;
        end
        if ~any(postMask)
            statusString = 'post_freq_not_found';
            return;
        end
        Gpost = mean(x.grangerPost(:,:,postMask),3,'omitnan');
        Gpost = removeSelfConnections(Gpost);
        outflowPost = mean(Gpost,2,'omitnan');
        inflowPost = mean(Gpost,1,'omitnan')';
        outflow = outflowPost;
        inflow = inflowPost;
    case {'diff','post-pre'}
        if ~preStatusOk || ~postStatusOk
            statusString = 'invalid_meta';
            return;
        end
        if ~any(preMask) || ~any(postMask)
            statusString = 'freq_not_found';
            return;
        end
        Gpre = mean(x.grangerPre(:,:,preMask),3,'omitnan');
        Gpost = mean(x.grangerPost(:,:,postMask),3,'omitnan');
        Gpre = removeSelfConnections(Gpre);
        Gpost = removeSelfConnections(Gpost);
        outflowPre = mean(Gpre,2,'omitnan');
        inflowPre = mean(Gpre,1,'omitnan')';
        outflowPost = mean(Gpost,2,'omitnan');
        inflowPost = mean(Gpost,1,'omitnan')';
        outflow = outflowPost - outflowPre;
        inflow = inflowPost - inflowPre;
    otherwise
        error('Unknown conditionType: %s',conditionType);
end
end

function G = removeSelfConnections(G)
n = size(G,1);
G(1:n+1:end) = NaN;
end

function selectionList = resolveSelection(selection,subjectNameList,groupNameList)
selectionText = char(selection);

if ismember(selectionText,subjectNameList)
    selectionList = {selectionText};
    return;
end

switch lower(selectionText)
    case 'meditator'
        selectionList = subjectNameList(strcmp(groupNameList,'Meditator'));
    case 'control'
        selectionList = subjectNameList(strcmp(groupNameList,'Control'));
    case 'all'
        selectionList = subjectNameList;
    otherwise
        error('Selection must be a subject code, ''Meditator'', ''Control'', or ''All''.');
end
end

function topoValues = summarizeMeasure(outflowAll,inflowAll,measureType,useMedianFlag)
switch lower(measureType)
    case 'outflow'
        valueMatrix = outflowAll;
    case 'inflow'
        valueMatrix = inflowAll;
    case 'netflow'
        valueMatrix = outflowAll - inflowAll;
    otherwise
        error('measureType must be ''outflow'', ''inflow'', or ''netflow''.');
end

if useMedianFlag
    topoValues = median(valueMatrix,2,'omitnan');
else
    topoValues = mean(valueMatrix,2,'omitnan');
end
end

function titleString = makeTitle(selection,protocolName,freqRange,conditionType,summaryMethod,numSubjects,measureType)
titleString = sprintf('Granger %s | %s | %s | %g-%g Hz | %s | %s of %d subjects', ...
    capitalizeLabel(measureType),char(selection),protocolName,freqRange(1),freqRange(2),upper(conditionType),summaryMethod,numSubjects);
end

function safeTag = sanitizeFileTag(selection)
safeTag = regexprep(char(selection),'[^a-zA-Z0-9]+','_');
end

function out = capitalizeLabel(in)
in = char(in);
if isempty(in)
    out = in;
else
    out = [upper(in(1)) lower(in(2:end))];
end
end

function statusString = getStatusString(metaStruct)
statusString = 'missing_meta';
if isfield(metaStruct,'status')
    statusString = metaStruct.status;
end
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

function protocolList = getDefaultProtocolList()
protocolList = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
end

function projectRoot = getProjectRoot()
projectRoot = fileparts(fileparts(mfilename('fullpath')));
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

function [montageChanlocs,montageLabels] = getMontageChanlocs(projectRoot)
layoutFolder = resolveExistingFolder({ ...
    fullfile(projectRoot,'Montages','Layouts','actiCap64_UOL'), ...
    fullfile(projectRoot,'Dependancies','Montages','Layouts','actiCap64_UOL')}, ...
    'actiCap64_UOL montage layout folder');
x = load(fullfile(layoutFolder,'actiCap64_UOL.mat'));
y = load(fullfile(layoutFolder,'actiCap64_UOLLabels.mat'));
montageChanlocs = x.chanlocs;
montageLabels = y.montageLabels(:,2);
end

function addEEGLABTopoplot(projectRoot)
if exist('topoplot','file') ~= 2
    eeglabFolder = resolveExistingFolder({ ...
        fullfile(projectRoot,'fieldtrip-20260211','external','eeglab'), ...
        fullfile(projectRoot,'Dependancies','fieldtrip-20260211','external','eeglab'), ...
        fullfile(projectRoot,'eeglab')}, ...
        'EEGLAB topoplot folder');
    addpath(eeglabFolder);
end
end

function folderName = resolveExistingFolder(candidateFolders,folderDescription)
folderName = '';
for iFolder = 1:numel(candidateFolders)
    if exist(candidateFolders{iFolder},'dir')
        folderName = candidateFolders{iFolder};
        return;
    end
end

error('Could not find %s. Checked:%s%s',folderDescription,newline,strjoin(candidateFolders,newline));
end

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
