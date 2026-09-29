folderName = fullfile(pwd);
pppp = genpath(folderName);
addpath(pppp)

beep off

set(0,'DefaultFigureColormap',colormap(jet(64)));

set(groot,'defaultLegendAutoUpdate','off')
set(0, 'DefaultFigureToolBar', 'figure');
set(0, 'DefaultAxesCreateFcn', @(src, ~) set(axtoolbar(src, 'default'), 'Visible', 'on'));

pname = '/Users/ivankova/Library/CloudStorage/GoogleDrive-irena.vanek@gmail.com/My Drive/Research/DOVuFRIS/Fris_Apres/code';

folderName = fullfile([pname '/tseries/sites/other_data/Fimbul/fimbul_front/FIMF/Bandwidth_subsets']);
pppp = genpath(folderName);
addpath(pppp)

folderName = fullfile([pname '/tseries/sites/other_data/Fimbul/fimbul_front/FIMF/FIMFY25']);
pppp = genpath(folderName);
addpath(pppp)

folderName = fullfile([pname '/Users/ivankova/Code/apres_ens/apres_ens/matlab/shared']);
pppp = genpath(folderName);
addpath(pppp)

