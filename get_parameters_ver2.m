function p = get_parameters_ver2()

% This is not quite same as in write_params_in, I changed just a few things. 
% Also I've used terser names here because
% there's going to be more typing on long lines.
pl = {'Restart from last run?'; 
    'End time of integration';
    'Time step';
    'Graphics averaging time';
    'Frequency of graphics output';
    'Interactive radiation?';
    'Interactive clouds?';
    'Interactive surface temp?';
    'Frequency of radiation calls';
    'Solar constant';
    'Latitude';
    'Starting month';
    'Starting day';
    'Starting hour';
    'Time-dependent radiation?';
    'Date-dependent radiation?';
    'Diurnal-average radiation?';
    'Annual-average radiation?';
    'Calculate ocean albedo?';
    'Surface albedo (if non-interactive)';
    'Mass concentration of CO2';
    'Mass concentration of CH4';
    'Mass concentration of N2O';
    'Mass concentration of CFC11';
    'Mass concentration of CFC12';
    'Interactive water vapor?';
    'Radiation H2O Multiplier';
    'Radiation O3 Multiplier';
    'Dry Adiabatic Adjustment';
    'Moist Convection';
    'Turbulent Fluxes';
    'Fraction of Surface Covered by Water';
    'Mixed layer depth';
    'Surface wind speed';
    'Cubic profile of omega?';
    'Extreme value of omega';
    'Period of omega';
    'P at which omega=0 (top)';
    'P at which omega=0 (bottom)';
    'P at which omega=extr. value';
    'Apply weak-temperature-gradient approximation?';
    'P above which sndg. fixed'};

% Also duplicated from write_params_in but modified slightly.
units = cell(size(pl));
units{2} = 'days';
units{3} = 'minutes';
units{4} = 'days';
units{5} = 'hours';
units{9} = 'hours';
units{10} = 'watts m^-2';
units{11} = 'degrees';
units{21} = 'ppm';
units{22} = 'ppm';
units{23} = 'ppb';
units{24} = 'ppt';
units{25} = 'ppt';
units{33} = 'm';
units{34} = 'm s^-1';
units{36} = 'mb hour^-1';
units{37} = 'days';
units{38} = 'mb';
units{39} = 'mb';
units{40} = 'mb';
units{42} = 'mb';

% cell of the default values.
p = {'n';
    100;
    5;
    10;
    12;
    'y';
    'y';
    'y';
    3;
    1360;
    15;
    3;
    1;
    0;
    'n';
    'n';
    'n';
    'y';
    'n';
    0.32;
    360;
    1.72;
    310;
    280;
    484;
    'y';
    1;
    1;
    'y';
    'y';
    'y';
    1;
    2;
    5; % wind speed
    % 10; % HW b)
    'y';
    0;
    2000;
    80;
    1000;
    750;
    'n';
    850};
%
j = 1;
while j == 1;
    menucell = cell(size(pl));
    for i = 1:length(pl) %Strictly the wrong way to do this.
        if(isempty(units{i}))
            menucell{i} = [num2str(pl{i}),' [',num2str(p{i}),']'];
        else
            menucell{i} = [num2str(pl{i}),' [',num2str(p{i}),' ',num2str(units{i}),']'];
        end
    end
    menucell{i+1} = 'Open Users Guide';
    menucell{i+2} = 'Run Model With Current Configuration';
    
    it = menu('Configure Model',menucell);

    if it == length(menucell)
        j = 0;
    elseif it ==  length(menucell)-1
        if isunix && ~ismac
             unix('pdfopen --file Users_guide.pdf');
        elseif ismac
             unix('open Users_guide.pdf');
        elseif ispc
             dos('start Users_guide.pdf');
        end     
    else
        if ischar(p{it})
            if strcmp(p{it},'y')
                p{it} = 'n';
            else
                p{it} = 'y';
            end
        else
            temp = inputdlg_new('New Value?');
            % First is to check if someone hit cancel, second is to check
            % if they hit OK without entering a value
            if ~isempty(temp) && ~isempty(temp{:})
                p{it} = str2num(cell2mat(temp));
            end
            
        end
    end
    
end













return
%{
% This maps the free parameters to the list of all the parameters.
% Parameters not on this list are not editable.

tabstop = 38; % number of spaces to align the default values to.

j = 1;
while j == 1
    clc
    fprintf(1,' Model Setup:\n\n');
    
    % Automatically generates menu text from the title of the parameters
    % and the current value. I realize that this is terrifying and
    % apologize.
    menutext = '';
    for i = 1:length(fp)
        % Terse names. We're going over just the free parameters,
        % displaying the title of the parameter pl{fp(i)}, and the value
        % it's currently set to p{fp(i)}.
        
        % That repmat stuff pads spaces so that the default values line up
        % in a nice column.
        
        % This is so the numbers are right-aligned.
        if i < 10
            numbertext = ['    ',num2str(i),') '];
        else
            numbertext = ['   ',num2str(i),') '];
        end
        if isempty(units{fp(i)})
            unitstring = '';
        else
            unitstring = [' ',num2str(units{fp(i)})];
        end
        menutext = [menutext,numbertext,pl{fp(i)},repmat(' ',1,tabstop-length(pl{fp(i)})-length(numbertext)),'[',num2str(p{fp(i)}),unitstring,']\n'];
    end
    menutext = [menutext,'\n    0) Run model with current configuration.\n\n'];
    
    fprintf(1,menutext);
    it = input('Enter value to change? ');
    if it==0 || it > length(fp)
        j = 0;
    else
        if ischar(p{fp(it)})
            if strcmp(p{fp(it)},'y')
                p{fp(it)} = 'n';
            else
                p{fp(it)} = 'y';
            end
        else
            p{fp(it)} = input('New Value? ');
            
        end
            
            
    end

end
%}