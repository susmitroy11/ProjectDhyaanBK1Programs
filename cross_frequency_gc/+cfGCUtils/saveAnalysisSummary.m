function saveAnalysisSummary(summaryFileName,analysisResult)
% saveAnalysisSummary Save an analysisResult struct after ensuring the folder exists.
analysisFolder = fileparts(summaryFileName);
cfGCUtils.ensureFolder(analysisFolder);
save(summaryFileName,'analysisResult','-v7.3');
end
