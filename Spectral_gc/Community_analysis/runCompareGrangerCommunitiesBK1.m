% runCompareGrangerCommunitiesBK1
% Convenience script to compare Granger communities between meditators and controls.

clear;
clc;

protocolList = {'M1','M2','EO1','EO2','EC1','EC2','G1','G2'};
badEyeCondition = 'ep';
badTrialVersion = 'v8';
saveOutputFlag = 1;
numPermutations = 5000;
conditionType = 'post'; % 'pre', 'post', or 'diff'
densityThreshold = 0.20;
gammaValue = 1;
numCommunityRuns = 50;

% Default bands used by the main comparison function:
% alpha     = 7-10 Hz
% beta      = 20-32 Hz
% highgamma = 30-80 Hz

bandInfo(1).name = 'theta-beta';
bandInfo(1).range = [4 16];

bandInfo(2).name = 'lowgamma';
bandInfo(2).range = [32 44];

bandInfo(3).name = 'highgamma';
bandInfo(3).range = [56 80];


statsUnpaired = cell(numel(protocolList),1);
statsPaired = cell(numel(protocolList),1);

for iProtocol = 1:numel(protocolList)
    protocolName = protocolList{iProtocol};
    fprintf('\nRunning community comparison for %s...\n',protocolName);

    statsUnpaired{iProtocol} = compareGrangerCommunitiesMeditatorsControlsBK1( ...
        protocolName,badEyeCondition,badTrialVersion,0,saveOutputFlag, ...
        numPermutations,bandInfo,conditionType,densityThreshold,gammaValue,numCommunityRuns);

    statsPaired{iProtocol} = compareGrangerCommunitiesMeditatorsControlsBK1( ...
        protocolName,badEyeCondition,badTrialVersion,1,saveOutputFlag, ...
        numPermutations,bandInfo,conditionType,densityThreshold,gammaValue,numCommunityRuns);
end

disp('Finished unpaired and paired Granger community comparisons for all 8 protocols.');
