%% Add LTVTools Source to MATLAB path
format compact;

% Add utils to path where roots are defined
addpath('./roots/');

% Add path to roots
addpath(tbxroot);
addpath(fullfile(tbxroot,'util'));
addpath(genpath(testroot));
addpath(genpath(demoroot));
addpath(helproot);