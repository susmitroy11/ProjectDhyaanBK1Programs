function analysisResult = runAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1()
% runAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1
%
% Interactive launcher for the main CF-GC workflows.
% Edit the configuration block below, then run this file from MATLAB to
% call one of the visualization or comparison entry points.
%
% Structure:
%   Option A: Visualization
%       Suboption 1: Single subject
%       Suboption 2: Multiple subjects
%   Option B: Comparison
%       Suboption 1: Single pair
%       Suboption 2: Multiple pairs

thisFolder = fileparts(mfilename('fullpath'));
analysisFolder = fullfile(thisFolder,'analysis');

%% Main selection
mainOption = 'comparison';
% 'visualization' or 'comparison'

subOption = 'multiple_pairs';
% For visualization:
%   'single_subject'
%   'multiple_subjects'
%
% For comparison:
%   'single_pair'
%   'multiple_pairs'

protocolNameList = 'all';
% 'all' or e.g. {'M1','M2'}

badEyeCondition = 'ep';
badTrialVersion = 'v8';
saveOutputFlag = 1;
forceRebuild = 0;

%% Option A inputs
singleSubjectInput = '019CKa';
multipleSubjectsInput = 'all';
% 'all' or e.g. {'013AR','064PK','019CKa'}

visualizationCfg = struct;
visualizationCfg.sampleCountThreshold = 20;
visualizationCfg.showSampleCountFigure = 1;
visualizationCfg.showOverviewFigure = 1;

%% Option B inputs
singlePairInput = {'013AR','064PK'};
% Can also be a pair index, e.g. 21

multiplePairsInput = 'all';
% 'all'
% or a vector of pair indices, e.g. [1 21]
% or a list of explicit pairs, e.g.
%   {{'013AR','064PK'}; {'019CKa','022SSP'}}

comparisonCfg = struct;
comparisonCfg.numPermutations = 5000;
comparisonCfg.minObservationCount = 20;
comparisonCfg.minValidSubjectsPerCell = 3;
comparisonCfg.fdrAlpha = 0.05;
comparisonCfg.makeOverviewFigure = 1;

%% Shared config
visualizationCfg.analysisFolder = analysisFolder;
visualizationCfg.saveOutputFlag = saveOutputFlag;
visualizationCfg.forceRebuild = forceRebuild;

comparisonCfg.analysisFolder = analysisFolder;
comparisonCfg.saveOutputFlag = saveOutputFlag;
comparisonCfg.forceRebuild = forceRebuild;

%% Run selected analysis
switch lower(mainOption)
    case 'visualization'
        switch lower(subOption)
            case 'single_subject'
                analysisResult = visualizeCrossFrequencyPhaseAmplitudeGrangerSingleSubjectBK1( ...
                    singleSubjectInput,protocolNameList,badEyeCondition,badTrialVersion,visualizationCfg);

            case 'multiple_subjects'
                analysisResult = visualizeCrossFrequencyPhaseAmplitudeGrangerMultipleSubjectsBK1( ...
                    multipleSubjectsInput,protocolNameList,badEyeCondition,badTrialVersion,visualizationCfg);

            otherwise
                error('For mainOption = visualization, subOption must be single_subject or multiple_subjects.');
        end

    case 'comparison'
        switch lower(subOption)
            case 'single_pair'
                analysisResult = compareCrossFrequencyPhaseAmplitudeGrangerSinglePairBK1( ...
                    singlePairInput,protocolNameList,badEyeCondition,badTrialVersion,comparisonCfg);

            case 'multiple_pairs'
                analysisResult = compareCrossFrequencyPhaseAmplitudeGrangerMultiplePairsBK1( ...
                    multiplePairsInput,protocolNameList,badEyeCondition,badTrialVersion,comparisonCfg);

            otherwise
                error('For mainOption = comparison, subOption must be single_pair or multiple_pairs.');
        end

    otherwise
        error('mainOption must be visualization or comparison.');
end

fprintf('Analysis finished. Outputs are organized under:\n%s\n',analysisFolder);
end
