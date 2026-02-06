function p = get_parameters()

% This is not quite same as in write_params_in, I changed just a few things. 
% Also I've used terser names here because
% there's going to be more typing on long lines.
pl = {'Restart from last run?';
    'Interactive radiation?';
    'Interactive clouds?';
    'Interactive surface temp?';
    'Solar constant';
    'Latitude';
    'Month';
    'Starting day';
    'Starting hour';
    'Time-dependent radiation?';
    'Date-dependent radiation?';
    'Diurnal-average radiation?';
    'Annual-average radiation?';
    'Calculate ocean albedo?';
    'Surface albedo';
    'Amount of CO2';
    'Mixed layer depth';
    'Surface wind speed';
    'Cubic profile of omega?';
    'Period of omega';
    'Extreme value of omega';
    'P at which omega=0';
    'P at which omega=extr. value';
    'P above which sndg. fixed';
    'Time step';
    'End time of integration';
    'Graphics averaging time';
    'Frequency of radiation calls';
    'Frequency of graphics output'};

% Also duplicated from write_params_in but modified slightly.
units = cell(29,1);
units{5} = 'watts/m^2';
units{6} = 'degrees';
units{8} = 'days';
units{9} = '1-24';
units{16} = 'ppm';
units{17} = 'm';
units{18} = 'm/s';
units{20} = 'days';
units{21} = 'mb/hour';
units{22} = 'mb';
units{23} = 'mb';
units{24} = 'mb';
units{25} = 'minutes';
units{26} = 'days';
units{27} = 'days';
units{28} = 'hours';
units{29} = 'hours';

p = {'n';
    'y';
    'y';
    'n';
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
    2;
    5;
    'y';
    2000;
    0;
    80;
    750;
    0;
    5;
    100;
    1;
    3;
    2};

% This maps the free parameters to the list of all the parameters.
% Parameters not on this list are not editable.
fp = [1:4,10:13,15,16,26:29];
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
            p{fp(it)} = inputdlg('New Value?');
            
        end
            
            
    end

end