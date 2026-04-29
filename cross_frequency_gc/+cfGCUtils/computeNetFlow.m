function netFlow = computeNetFlow(matrixData)
% computeNetFlow Convert a directed connectivity matrix into net flow.
outFlow = mean(matrixData,2,'omitnan');
inFlow = mean(matrixData,1,'omitnan')';
netFlow = outFlow - inFlow;
end
