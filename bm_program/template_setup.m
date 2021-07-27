function [setup_id] =  template_setup(setup)
%SETUP FILE
disp(setup)
% =========================================================================
%parameter setup: names of varied paramter have to be defined below with range
IDENTIFIER = 'test';

varied_parameter{1} = 'ecen';
varied_parameter{2} = 'lin_system_tol';

parameter_combination = 'mesh';


% =========================================================================
% setup parameter of setup - those typical for the used program (CIM RIM)




parameter = vars2struct();


% =========================================================================
% HPC run_preparation (general)
% GET range
parameter_range = cell(numel(varied_parameter),1);
for k = 1:numel(varied_parameter)
    parameter_range{k,1} = eval(varied_parameter{k});
end

% GET range_cell (list of actual condor / slurm runs)
switch parameter_combination 
    case 'mesh'
        switch numel(varied_parameter)
            case 1
                if iscell(parameter_range{1})
                    parameter_vector_4_run = parameter_range{1}(:);
                else
                    parameter_vector_4_run = num2cell(parameter_range{2}(:));
                end
            case 2
                [xx{1},xx{2}] = meshgrid(parameter_range{1}, parameter_range{2});
                parameter_vector_4_run = cell(numel(xx{1}), numel(varied_parameter));
                for k = 1:2
                    if iscell(xx{k})
                        parameter_vector_4_run(:,k) = xx{k}(:);
                    else
                        parameter_vector_4_run(:,k) = num2cell(xx{k}(:));
                    end
                end
               
            case 3
                [xx{1},xx{2},xx{3}] = meshgrid(parameter_range{1},parameter_range{2},parameter_range{3});
                parameter_vector_4_run = cell(numel(xx{1}), numel(varied_parameter));
                for k = 1:3
                    if iscell(xx{k})
                        parameter_vector_4_run(:,k) = xx{k}(:);
                    else
                        parameter_vector_4_run(:,k) = num2cell(xx{k}(:));
                    end
                end

            otherwise
                error('Mesh functionality for more than 3 dimensions not yet implemented')
        end
        
        
    case 'serial'
        if any(cellfun(@(x) numel(x), parameter_range) ~= numel(parameter_range{1}))
            error('For serial combinations of parameters the number of variations for each parameter have to be the same')
        end
        parameter_vector_4_run = cell(numel(parameter_range{1}), numel(varied_parameter));
        for k = 1:numel(varied_parameter)
            if iscell(parameter_range{k}(:))
                parameter_vector_4_run(:,k) = parameter_range{k}(:);
            else
                parameter_vector_4_run(:,k) = num2cell(parameter_range{k});
            end
        end
    case 'fractional'
        %not yet implemented "Versuchsplan", design of experiments
        %factorial experiment
end


date_format = 'yymmdd';
%time_format = 'HH:MM';  % optional
date_string = datestr(now, [date_format]);

%check if setup_file already exists:
listing = dir(['./data/setup_', IDENTIFIER, '_', date_string, '_*.mat']);
number = numel(listing);

setup_id = [IDENTIFIER, '_', date_string, '_', num2str(number)];

data_directory = ['./data/data_', setup_id];

par_setup = vars2struct('IDENTIFIER', 'setup_id', 'data_directory', 'parameter_combination', 'parameter_vector_4_run', 'parameter_range', 'varied_parameter');
% =========================================================================

% ------------------------------ save -------------------------------------


setup = vars2struct('parameter', 'par_setup');

path_to_save = './data/';
save([path_to_save, 'setup_', setup_id ,'.mat'], 'setup')

% %check compile date
% if ~isdeployed()
%     GROUNDPATH = which('contour.m');
%     change_list = dir(GROUNDPATH);
%     run_folder = fileparts(fileparts(mfilename('fullpath'))); 
%     compile_list = dir([run_folder, '/contour_program/contour']);
%     if isempty(compile_list)
%         msgbox('Warning: Program not compiled yet.')
%     else
%         if datetime(change_list.datenum,  'ConvertFrom','datenum') > datetime(compile_list.datenum, 'ConvertFrom','datenum')
%             msgbox('Warning: There were changes after compiling', 'Compiled?', 'warn')
%         end
%     end
%     disp(['Compile date: ', compile_list.date])
% end
disp(['Created file name: ', path_to_save, 'setup_', setup_id ,'.mat', ', number of parameter runs: ', num2str(size(parameter_vector_4_run,1))])

end