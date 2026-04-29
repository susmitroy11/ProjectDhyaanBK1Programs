function ensureFolder(folderName)
% ensureFolder Create a folder if it does not already exist.
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
