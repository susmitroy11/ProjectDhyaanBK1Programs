function [gridData,gridSampleCount,pairInfo,phaseBandList,ampBandList] = unpackBandPairGrid(results)
% unpackBandPairGrid Extract the requested phase-amplitude subgrid from results.
numROIs = numel(results.roiLabels);
phaseBandList = unique(results.targetBandPairs(:,1),'stable');
ampBandList = unique(results.targetBandPairs(:,2),'stable');

numPhaseBands = numel(phaseBandList);
numAmpBands = numel(ampBandList);

gridData = nan(numROIs,numROIs,numPhaseBands,numAmpBands);
gridSampleCount = nan(numROIs,numROIs,numPhaseBands,numAmpBands);
pairInfo = repmat(struct('phaseBandIndex',[],'ampBandIndex',[], ...
    'phaseRange',[],'ampRange',[],'label',''),numPhaseBands,numAmpBands);

for iPhase = 1:numPhaseBands
    for iAmp = 1:numAmpBands
        phaseBandIndex = phaseBandList(iPhase);
        ampBandIndex = ampBandList(iAmp);

        gridData(:,:,iPhase,iAmp) = results.phaseToAmpGC(:,:,phaseBandIndex,ampBandIndex);
        gridSampleCount(:,:,iPhase,iAmp) = results.sampleCount(:,:,phaseBandIndex,ampBandIndex);

        phaseRange = results.bandInfo(phaseBandIndex).range;
        ampRange = results.bandInfo(ampBandIndex).range;
        pairInfo(iPhase,iAmp).phaseBandIndex = phaseBandIndex;
        pairInfo(iPhase,iAmp).ampBandIndex = ampBandIndex;
        pairInfo(iPhase,iAmp).phaseRange = phaseRange;
        pairInfo(iPhase,iAmp).ampRange = ampRange;
        pairInfo(iPhase,iAmp).label = sprintf('Phase %g-%g Hz -> Amp %g-%g Hz', ...
            phaseRange(1),phaseRange(2),ampRange(1),ampRange(2));
    end
end
end
