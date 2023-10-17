% COMSOL with MATLAB - Simulink
% Adjust active recording electrode and grey matter conductivity value
% Generates COMSOL files with solutions
% Lauren R. Madden | 2023

% import Simulink
import com.comsol.model.*
import com.comsol.model.util.*

% set up directories
cd ..\COMSOL_files;

fileprefix = 'Feline_S12_DRG';

% load basic model
model = mphload('Feline_S12_DRG_ppd');

elec_range = 36:2:46; % electrode contact number range
elec_boundary = 77:-2:67; % COMSOL's boundary identifcation numbers for each electrode
sigma_range = [0.083,0.26,0.57]; % range of conductivity values for grey matter (Ranck and BeMent 1965)

% show progress bar
ModelUtil.showProgress(true);

% per electrode contact
for e = 1:size(elec_range,2)
    % per sigma value
    for s = 1:size(sigma_range,2)

        elec = elec_range(e);
        elec_str = string(elec);
        sigma = sigma_range(s);
        sigma_str = strrep(string(sigma),'.','p');
        
        % set grey matter conductivity & resistivity
        model.component('comp1').material('mat7').propertyGroup('def').set('electricconductivity',sigma);
        model.component('comp1').material('mat7').propertyGroup('def').set('resistivity',1/sigma);
        
        if s == 1 % if first sigma value for new electrode contact loop, change which electrode contact is active and finely meshed
            
            % disable active contact's Floating Potential
            model.component('comp1').physics('ec').feature(strcat('fp',elec_str)).active(false);
            % set Terminal to active contact
            model.component('comp1').physics('ec').feature('term1').selection.set([elec_boundary(e)]);
            % enable previous contact's Floating Potential
            if e ~= 1 % if contact is not first one
                model.component('comp1').physics('ec').feature(strcat('fp',string(elec_range(e-1)))).active(true);
            end
    
            % set meshing concentration to active electrode contact's Boundary
            model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').selection.geom('geom1', 2);
            model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').selection.set([elec_boundary(e)]);
            
            % mesh
            model.component('comp1').mesh('mesh1').run;
            disp(strcat(fileprefix,'_e',elec_str,'_sigma',sigma_str,' is meshed...'))

        end
        
        % solve
        model.sol('sol1').runAll;
        disp(strcat(fileprefix,'_e',elec_str,'_sigma',sigma_str,' is solved...'))
        
        % save solution
        mphsave(model,strcat(fileprefix,'_e',elec_str,'_sigma',sigma_str))

        disp(strcat(fileprefix,'_e',elec_str,'_sigma',sigma_str,' is saved.'))

    end
end
