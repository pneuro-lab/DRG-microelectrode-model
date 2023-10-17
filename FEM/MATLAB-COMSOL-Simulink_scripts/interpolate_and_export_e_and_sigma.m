% COMSOL with MATLAB - Simulink
% Adjust active recording electrode and grey matter conductivity value
% Interpolates NEURON compartment center coordinates with COMSOL solutions
% Lauren R. Madden | 2023

% ------ edit this for desired NEURON model center coordinate file-------
nrn_center_file = 'nrn_centers.txt';
%------------------------------------------------------------------------

% import Simulink
import com.comsol.model.*
import com.comsol.model.util.*

% set up directories
cd ..\..\
cur_dir = pwd;
model_dir = 'FEM\COMSOL_files\';
nrn_center_files = '\InterpolationFiles\';
interp_save_dir = '\InterpolationFiles\';


elec_range = 36:2:46; % electrode contact number range
elec_boundary = 77:-2:67; % COMSOL's boundary identifcation numbers for each electrode
sigma_range = ["083","26","57"]; % range of conductivity values for grey matter (Ranck and BeMent 1965)

fileprefix = 'Feline_S12_DRG';

% per electrode contact
for e = 1:size(elec_range,2)

    elec = elec_range(e);
    elec_str = string(elec);
    
    % per sigma value
    for s = 1:size(sigma_range,2)
        sigma_str = sigma_range(s);
        
        ncfsuffix = split(nrn_center_file,'.');
        ncfsuffix = split(ncfsuffix{1},'_');
        if ~strcmp(ncfsuffix{end},'centers')
            ncfsuffix = strcat('_',ncfsuffix{end});
        else
            ncfsuffix = "";
        end

        % load model with solution
        model = mphload(strcat(model_dir,fileprefix,'_e',elec_str,'_sigma0p',sigma_str));
        
        model.result.export('data1').setIndex('expr', 'V', 0);
        model.result.export('data1').set('filename',strcat(cur_dir,interp_save_dir,'interpolation_e',elec_str,'_sigma0p',sigma_str,ncfsuffix,'.txt'));
        
        % set export data format with struct property
        model.result.export('data1').set('coordfilename', strcat(cur_dir,nrn_center_files,nrn_center_file));
        
        % run export
        model.result.export('data1').run;
        disp(strcat('interpolation_e',elec_str,'_sigma0p',sigma_str,ncfsuffix,' is saved'));
        
        
    end
end
