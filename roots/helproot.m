function out = helproot()
% Returns Help Root
out = fullfile(fileparts(fileparts(mfilename('fullpath'))),'doc');