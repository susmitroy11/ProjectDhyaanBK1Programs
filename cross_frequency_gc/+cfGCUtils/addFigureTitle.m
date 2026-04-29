function addFigureTitle(hFig,titleText)
% addFigureTitle Add a figure-wide title with backward compatibility.
if exist('sgtitle','file') == 2
    sgtitle(titleText);
else
    annotation(hFig,'textbox',[0 0.965 1 0.03],'String',titleText, ...
        'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');
end
end
