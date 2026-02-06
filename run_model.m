% run_model.m
% Notes:
% - The Fortran executable `rc_ver2*` reads control parameters from
%   `params_ver2.in` and initial conditions from `sounding.in`.
% - This wrapper uses `get_parameters_ver2()` + `write_params_in_ver2()` to
%   generate `params_ver2.in`, and `write_sounding_in()` to generate
%   `sounding.in` (including SST).
%                Initial souding is generated in this code -> write_sounding_in.m (SST setting)

% Runs the model and plots the output. 
clear
clc
%  Allows used to specify output directory name. (Output is written to
%  'output' but directory is then renamed)
%
%dirname=cell2mat(inputdlg_new('Enter directory name (Default="output1")'));
%if isempty(dirname)
    dirname='output1';
%end    
%unix('rm ./output/*');
parameters = get_parameters_ver2();

% Are we resuming a run? If so we need to copy sounding.out
%if strcmp(parameters{1},'y') 
%    % Does the file exist? I actually check for both the output dir and the
%    if (exist(dirname,'dir')==7) && (exist([dirname,'/sounding.out'],'file')==2)
%        copyfile([dirname,'/sounding.out'],'./sounding.out');
%    else
%        error('No old run in the output directory to resume!')
%    end
%end

write_params_in_ver2(parameters);


% Change sounding if needed (SST, etc.)
write_sounding_in();

mkdir('output');
if isunix && ~ismac
    %cmd='mkdir -p output';
    %[status,result]=unix(cmd);    
%    unix('chmod +x rc_ver2unix');
     unix('./rc_ver2unix');
%    unix('./a.out');
elseif ismac
    if exist('rc_ver2mac','file')==2
        unix('chmod +x rc_ver2mac');
        unix('./rc_ver2mac');
    else
        error('Mac version must be compiled first! Use gfortran -o rc_ver2mac rc_ver2.f');
    end
elseif ispc
    %cmd='mkdir output';
    %[status,result]=dos(cmd);
    dos('rc_ver2dos.exe');
end
rcm; % plot scripts

% This is more platform-independent.
if ~strcmp(dirname,'output')
    a = rmdir(dirname,'s');
    a = movefile('output',dirname);
end
%{
if isunix && ~ismac && ~strcmp(dirname,'output')
    filerm=['rm -f -r ' dirname];
    [status,result]=unix(filerm);
    filecmd=['mv output ' dirname];
    unix(filecmd);
elseif ispc && ~strcmp(dirname,'output')
    filerm=['rmdir /S/Q ' dirname];
    [status,result]=dos(filerm);
    filecmd=['rename output ' dirname];
    dos(filecmd);
end
%}