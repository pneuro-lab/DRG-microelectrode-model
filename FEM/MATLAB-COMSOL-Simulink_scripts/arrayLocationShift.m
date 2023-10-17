% COMSOL with MATLAB - Simulink
% Shift electrode array position in X or Z direction
% Generates COMSOL files with solutions
% Lauren R. Madden | 2023

%---------------------- comment out one option ---------------------------
shift_type = 'Xshift'; % shift array location along the axis of peripheral/central branches
% shift_type = 'Zshift'; % shift array location perpendicular to the axis of the peripheral/central branches

%---------------------- comment out one option ---------------------------
shift_subtype = 'ppd'; % contacts of array face paralell to axis of peripheral/central branches
% shift_subtype = 'pll'; % contacts of array face the central branch
%-------------------------------------------------------------------------

% import Simulink
import com.comsol.model.*
import com.comsol.model.util.*

% set up directories
cd ../COMSOL_files;

elec_range = 36:2:46; % electrode contact number range
if strcmp(shift_subtype,'ppd')
    elec_boundary = 77:-2:67; % COMSOL's boundary identifcation numbers for each electrode - ppd
else
    elec_boundary = 103:-1:98; % COMSOL's boundary identifcation numbers for each electrode - pll
end

% set shift values depending on axis direction
if strcmp('Xshift',shift_type)
    shift_range = [0, 0.01, 0.035, 0.085, 0.185, 0.285]; % array lead body shift values in mm, in positive X direction (along central axon)
else
    shift_range = [0,-0.01,-0.035,-0.085,-0.185,-0.280]; % array lead body shift values in mm, in negative Z direction (shifting laterally), Distance greater than 0.281 mm causes geometry error
end

% set Grey Matter conductivity & resistivity
model.component('comp1').material('mat7').propertyGroup('def').set('electricconductivity',0.26);
model.component('comp1').material('mat7').propertyGroup('def').set('resistivity',1/0.26);

% show progress bar
ModelUtil.showProgress(true);

% per electrode contact
for e = 1:size(elec_range,2)

    % load basic model
    model = mphload(['Feline_S12_DRG_',shift_subtype]);
    
    elec = elec_range(e);
    elec_str = string(elec);

    % disable active contact's Floating Potential
    model.component('comp1').physics('ec').feature(strcat('fp',elec_str)).active(false);
    % set Terminal to active contact
    model.component('comp1').physics('ec').feature('term1').selection.set([elec_boundary(e)]);
    % enable previous contact's Floating Potential if contact is not first one
    if e ~= 1
        model.component('comp1').physics('ec').feature(strcat('fp',string(elec_range(e-1)))).active(true);
    end

    % set meshing concentration to active electrode contact's Boundary
    model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').selection.geom('geom1', 2);
    model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').selection.set([elec_boundary(e)]);
    
    % per lead position shift value
    for s = 1:size(shift_range,2)

        shift = shift_range(s);
        shift_str = strrep(strrep(string(shift),'.','p'),'-','n');

        % set shift parameter to correct value
        if strcmp('Xshift',shift_type)
            model.param.set('shiftx',shift);
            model.param.set('shiftz',0);
        else
            model.param.set('shiftz',shift);
            model.param.set('shiftx',0);
        end
        
        % mesh
        model.component('comp1').mesh('mesh1').run;
        disp(strcat('Feline_S12_DRG_',shift_subtype,'_e',elec_str,'_',shift_type,shift_str,' is meshed...'))
        
        % molve
        model.sol('sol1').runAll;
        disp(strcat('Feline_S12_DRG_',shift_subtype,'_e',elec_str,'_',shift_type,shift_str,' is solved...'))
        
        % save solution
        mphsave(model,strcat('Feline_S12_DRG_',shift_subtype,'_e',elec_str,'_',shift_type,shift_str))

        disp(strcat('Feline_S12_DRG_',shift_subtype,'_e',elec_str,'_',shift_type,shift_str,' is saved.'))

    end
end

