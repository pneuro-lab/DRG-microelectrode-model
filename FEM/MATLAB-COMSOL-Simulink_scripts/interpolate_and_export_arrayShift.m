% COMSOL with MATLAB - Simulink
% Interpolate and Export: Shift electrode array position in X or Z direction
% Lauren R. Madden | 2023

%---- comment out one option
shift_type = 'Xshift'; % shift array location along the axis of peripheral/central branches
% shift_type = 'Zshift'; % shift array location perpendicular to the axis of the peripheral/central branches

%---- comment out one option
shift_subtype = 'ppd'; % contacts of array face paralell to axis of peripheral/central branches
% shift_subtype = 'pll'; % contacts of array face the central branch

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

% Set shift values depending on axis direction
if strcmp('Xshift',shift_type)
    shift_range = [0, 0.01, 0.035, 0.085, 0.185, 0.285]; % array lead body shift values in mm, in positive X direction (along central axon)
else
    shift_range = [0,-0.01,-0.035,-0.085,-0.185,-0.285]; % array lead body shift values in mm, in negative Z direction (shifting laterally)
end

for e = 1:size(elec_range,2) % per electrode
    
    elec = elec_range(e);
    elec_str = string(elec);
    
    for s = 1:size(shift_range,2) % per lead position shift value

        shift = shift_range(s);
        shift_str = strrep(strrep(string(shift),'.','p'),['' '-'],'n');

        ncfsuffix = split(nrn_center_file,'.');
        ncfsuffix = split(ncfsuffix{1},'_');
        if ~strcmp(ncfsuffix{end},'centers')
            ncfsuffix = strcat('_',ncfsuffix{end});
        else
            ncfsuffix = "";
        end

        % Open Model
        model = mphload(strcat(model_dir,'Feline_S12_DRG_',shift_subtype,'_e',elec_str,'_',shift_type,shift_str));
        
        % set export expression and filename
        model.result.export('data1').setIndex('expr', 'V', 0);
        model.result.export('data1').set('filename',strcat(cur_dir,interp_save_dir,'interpolation_',shift_subtype,'_e',elec_str,'_',shift_type,shift_str,ncfsuffix,'.txt'));
        
        % set export data format with struct property
        model.result.export('data1').set('coordfilename', strcat(cur_dir,nrn_center_files,nrn_center_file));

        % Run export
        model.result.export('data1').run;

        % progress printout
        disp(strcat('interpolation_',shift_subtype,'_e',elec_str,'_',shift_type,shift_str,ncfsuffix,'.txt is saved'));

    end
end


