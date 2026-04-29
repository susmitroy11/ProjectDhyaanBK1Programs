function runSaveGrangerData(saveFTDataFlag)
% runSaveGrangerData Batch runner for BK1 Granger causality analysis.
%
% It reuses the existing saveFTData pipeline and then saves one Granger file
% per subject and protocol under connectivityProjectCodes/savedDataGranger.

if ~exist('saveFTDataFlag','var') || isempty(saveFTDataFlag)
    saveFTDataFlag = 0;
end

folderSourceString = getProjectRoot();
ensureBK1Paths(folderSourceString);
ensureFieldTripOnPath(folderSourceString);

goodSubjectList = getGoodSubjectsBK1;
[allSubjectNames,expDateList] = getDemographicDetails('BK1');

ftDataFolder = fullfile(folderSourceString,'data','ftData');
ensureFolder(ftDataFolder);

badEyeCondition = 'ep';
badTrialVersion = 'v8';
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
useTheseIndices = 1:length(goodSubjectList);

if saveFTDataFlag
    for iSubject = 1:length(useTheseIndices)
        subjectName = goodSubjectList{useTheseIndices(iSubject)};
        disp(['Saving FT data for subject ' subjectName]);
        expDate = expDateList{strcmp(subjectName,allSubjectNames)};
        saveFTData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,ftDataFolder);
    end
end

stRange = [0.25 1.25];
gcParams = struct;
gcParams.modelOrder = 10;
gcParams.mvarToolbox = 'biosig';
gcParams.maxFreq = 100;
gcParams.freqStep = 1;
gcParams.feedback = 'none';
gcParams.savedDataFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataGranger');

for iSubject = 1:length(useTheseIndices)
    subjectName = goodSubjectList{useTheseIndices(iSubject)};
    disp(['Saving Granger data for subject ' subjectName]);
    saveGrangerData(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,stRange,gcParams);
end
end

function folderSourceString = getProjectRoot()
folderSourceString = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

function ensureBK1Paths(folderSourceString)
addpath(genpath(fullfile(folderSourceString,'ProjectDhyaanBK1Programs')));
addpath(genpath(fullfile(folderSourceString,'CommonPrograms')));
addpath(genpath(fullfile(folderSourceString,'Montages')));
end

function ensureFieldTripOnPath(folderSourceString)
fieldTripPath = fullfile(folderSourceString,'fieldtrip-20260211');
addpath(fieldTripPath);
ft_defaults;
end

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
