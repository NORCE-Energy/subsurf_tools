
% Path information

toolbox_dir = '/media/tubh/IntData/Bhakta/NIOR Project/Work in progress/MatlabCodes';
IFT_dir = [toolbox_dir '/rfmatlab/main']; % IRISFilterToolbox directory
work_dir = pwd; 

options.work_dir = work_dir;

cd(IFT_dir)
cd ..
tmp_dir = pwd;
warning off
rmpath(genpath(tmp_dir)); % this path information is used to clean up those in MATLAB search path to avoid miscalling overloaded scripts
warning on
needed_subfolder = {'eclipse','eclipseKalman','iofun','Kalman',...
                    'plot2d','postProcess','preProcess','stats','seismicHM/NORNE-2D/local_scripts'};
for i = 1 : length(needed_subfolder)
    addpath([IFT_dir '/' needed_subfolder{i}]);
end


cd(work_dir);
close all;
