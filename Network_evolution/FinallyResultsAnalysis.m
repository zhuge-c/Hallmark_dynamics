%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finaly results analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

directory = '../../';

%  dir: get all the files under directory
folders = dir(directory);

% screen the folder 
folderNames = {folders([folders.isdir]).name};

