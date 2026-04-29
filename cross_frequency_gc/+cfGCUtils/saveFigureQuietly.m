function saveFigureQuietly(hFig,fileName)
% saveFigureQuietly Save a figure and close it automatically afterward.
cfGCUtils.ensureFolder(fileparts(fileName));
cleanObj = onCleanup(@() closeFigureIfValid(hFig)); %#ok<NASGU>

set(hFig,'PaperPositionMode','auto');
drawnow;

[~,~,ext] = fileparts(fileName);
ext = lower(ext);

try
    if exist('exportgraphics','file') == 2
        switch ext
            case '.png'
                exportgraphics(hFig,fileName,'Resolution',200);
            case {'.pdf','.eps','.svg'}
                exportgraphics(hFig,fileName,'ContentType','vector');
            otherwise
                exportgraphics(hFig,fileName);
        end
        return;
    end
catch
end

try
    switch ext
        case '.png'
            print(hFig,fileName,'-dpng','-r200');
        case '.pdf'
            print(hFig,fileName,'-dpdf','-painters');
        case '.eps'
            print(hFig,fileName,'-depsc','-painters');
        otherwise
            saveas(hFig,fileName);
    end
    return;
catch
end

saveas(hFig,fileName);
end

function closeFigureIfValid(hFig)
if ishghandle(hFig)
    close(hFig);
end
end
