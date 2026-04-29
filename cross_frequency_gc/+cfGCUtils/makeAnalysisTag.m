function tag = makeAnalysisTag(baseLabel,optionsText)
% makeAnalysisTag Build a stable folder tag from a label and option hash.
hashText = cfGCUtils.makeOptionsHash(optionsText);
tag = [cfGCUtils.sanitizePathToken(baseLabel) '__' hashText(1:10)];
end
