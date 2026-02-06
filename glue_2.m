% glue.m

% Runs the model and plots the output. 
clear
clc

%unix('rm ./output/*');
parameters = get_parameters();

write_params_in(parameters);

% It never changes.
% write_sounding_in();
if isunix && ~ismac
    unix('./rcnewrad');
elseif ispc
    dos('rcnew.exe');
end
rcmenu_new;