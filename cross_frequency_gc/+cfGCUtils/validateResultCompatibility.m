function validateResultCompatibility(preResults,postResults)
% validateResultCompatibility Confirm that two CF-GC result structs align.
if ~isequal(preResults.roiLabels,postResults.roiLabels)
    error('PRE and POST ROI labels do not match.');
end

if ~isequal(preResults.targetBandPairs,postResults.targetBandPairs)
    error('PRE and POST targetBandPairs do not match.');
end

if ~isequal(size(preResults.phaseToAmpGC),size(postResults.phaseToAmpGC))
    error('PRE and POST phaseToAmpGC sizes do not match.');
end
end
