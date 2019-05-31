% Run this code when you extracted the "GC Model - full.zip" 
% including this script at its final loation.

p = fileparts(mfilename('fullpath')); % get path to this script
addpath(genpath(p)); % add folder including subfolders to the Matlab path
fprintf('Added %s and its subdirectories to the Matlab path!\n',p);
