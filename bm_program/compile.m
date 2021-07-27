
% find out which system we are working on (VSC, ITP, Laptop)
%NOTE: Now on all Unix systems same Groundpath!

if ismac
    % Code to run on Mac plaform
elseif isunix
    % Code to run on Linux plaform
    GROUNDPATH = '~/QUENTIN/bm_program/';
    % GROUNDPATH = '/temp/dorn/run/bm_program/';
elseif ispc
    % Code to run on Windows platform
    
    GROUNDPATH = 'C:/QUENTIN/bm_program/';
else
    disp('Platform not supported')
end


%Alternative
% comp = computer;
% if contains(comp,'LN') 
%     [~,hostname] = system('hostname');
%     if contains(hostname, 'faep')
%         GROUNDPATH = '/temp/dorn/run/bm_program/';
%     else %vsc
%         GROUNDPATH = '~/Matlab/run/bm_program/';
%     end
% elseif contains(comp, 'PC')
%     GROUNDPATH = 'C:/Matlab/run/bm_program/';
% end



% compile C - function
%cd ([GROUNDPATH, 'functions', SLASH, 'c_functions', SLASH])

%mex GCC='/usr/local/bin/gcc-4.7.4' COMPFLAGS='$COMPFLAGS -Wall' c_GF.c
%mex COMPFLAGS='$COMPFLAGS -Wall' c_GF.c
current_path = mfilename('fullpath');
cd(current_path(1:end-length(mfilename)))


% compile program
%GROUNDPATH = '~/bm_program/';
addpath(GROUNDPATH)
addpath([GROUNDPATH, 'classes/'])
addpath([GROUNDPATH, 'functions/'])
addpath([GROUNDPATH, 'functions/class_functions/'])
addpath([GROUNDPATH, 'functions/c_functions/'])
addpath([GROUNDPATH, 'functions/c_functions/m_GF/'])
addpath([GROUNDPATH, 'functions/plot/'])
    
mcc -m -R -nojvm exec_bm.m

