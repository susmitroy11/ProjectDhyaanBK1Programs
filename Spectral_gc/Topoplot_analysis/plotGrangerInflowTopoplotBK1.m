function topo = plotGrangerInflowTopoplotBK1(selection,protocolName,freqRange,conditionType,badEyeCondition,badTrialVersion,saveFigureFlag,useMedianFlag)
% plotGrangerInflowTopoplotBK1 Plot BK1 scalp topography of Granger inflow.
% This is a thin wrapper around plotGrangerOutflowTopoplotBK1 that switches
% the measure type to inflow.
%
% Example:
%   topo = plotGrangerInflowTopoplotBK1('Meditator','M1',[30 80],'diff','ep','v8',1);
%   topo = plotGrangerInflowTopoplotBK1('BothGroups','all',[30 80],'all','ep','v8',1);

if ~exist('selection','var');         selection = []; end
if ~exist('protocolName','var');      protocolName = []; end
if ~exist('freqRange','var');         freqRange = []; end
if ~exist('conditionType','var');     conditionType = []; end
if ~exist('badEyeCondition','var');   badEyeCondition = []; end
if ~exist('badTrialVersion','var');   badTrialVersion = []; end
if ~exist('saveFigureFlag','var');    saveFigureFlag = []; end
if ~exist('useMedianFlag','var');     useMedianFlag = []; end

topo = plotGrangerOutflowTopoplotBK1(selection,protocolName,freqRange,conditionType,badEyeCondition,badTrialVersion,saveFigureFlag,useMedianFlag,'inflow');
end
