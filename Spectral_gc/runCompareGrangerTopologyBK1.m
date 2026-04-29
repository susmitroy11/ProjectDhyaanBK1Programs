% runCompareGrangerTopologyBK1
% Convenience script to compare Granger topology between meditators and controls.

clear;
clc;

protocolList = {'M1','M2','EO1','EO2','EC1','EC2','G1','G2'};
badEyeCondition = 'ep';
badTrialVersion = 'v8';
saveOutputFlag = 1;
numPermutations = 5000;
conditionType = 'diff'; % 'pre', 'post', or 'diff' (post - pre)
hubZThreshold = 1;

% Default bands used by the main comparison function:
% alpha     = 7-10 Hz
% beta      = 20-32 Hz
% highgamma = 30-80 Hz
bandInfo = [];

statsUnpaired = cell(numel(protocolList),1);
statsPaired = cell(numel(protocolList),1);

for iProtocol = 1:numel(protocolList)
    protocolName = protocolList{iProtocol};
    fprintf('\nRunning topology comparison for %s...\n',protocolName);

    statsUnpaired{iProtocol} = compareGrangerTopologyMeditatorsControlsBK1( ...
        protocolName,badEyeCondition,badTrialVersion,0,saveOutputFlag, ...
        numPermutations,bandInfo,conditionType,hubZThreshold);

    statsPaired{iProtocol} = compareGrangerTopologyMeditatorsControlsBK1( ...
        protocolName,badEyeCondition,badTrialVersion,1,saveOutputFlag, ...
        numPermutations,bandInfo,conditionType,hubZThreshold);
end

disp('Finished unpaired and paired Granger topology comparisons for all 8 protocols.');
