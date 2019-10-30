% PhenoCAT project, run script PER PLATE

which cropNuclei_saveRG_batches
which IllumCorrect
which fixNonNumericalValueInImage

% CHANGE PLATE NUMBER! (CP611 - CP616)
pathToData = '/Volumes/MotzBook/2018_12_19_MOTCscreenData/MOTC_screen_01/20151103_CN_CP614-1aa';

% CHANGE FOLDER (plate1-plate6)
pathToResults =  '/Volumes/MotzBook/2018_12_17_deep_learning/plate4';

load('/Users/iMotz/PhenoCAT_private/MATLAB_image_preparation/metaData.mat');

n_batches = 10;
channels = {'A01Z01C03','A02Z01C02'}; % R - G ; have to include action numbers
clims = [100 1000;100 5000]; % R - G
perimeter = 2^8;

cropNuclei_saveRG_batches(pathToData,metaData,pathToResults,channels,clims,perimeter,n_batches)

