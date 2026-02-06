function write_params_in(defaults)

b = fopen('params.in','w');

% This is what we call programming by rote.
firstlines = '  Parameter                       Value          Units\n  ---------                       -----          -----\n';
paramlines = {'Restart from last run?';
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
    'Surface albedo (if ''n'' above)';
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

% This is going to get stuck in a different file
%{
defaults = {'n';
    'y';
    'y';
    'y';
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
    10;
    3;
    12};
%}
% The units aren't defined for everything.
units = cell(29,1);
units{5} = 'watts/m^2 '; % The space was there in the original. I suspect this isn't that sensitive to whitespace but I'm reproducing the template exactly.
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

% We need a way of showing that month and starting day are integers aligned
% with the decimal point, but all the other numbers aren't;
is_value_int = zeros(size(defaults));
is_value_int(7:8) = 1;

% So if we're writing a float, we pad with spaces up to column 33. If
% we're writing a bool (y/n) or an integer, we'll need to add four more
% spaces to bring it in line with where the decimal point is in the float.
tabstop = 32; 

% I will be referring to this bit of code below. It should be a function
% but I wanted to remain a single file. To produce arbitrary numbers of
% spaces:
% repmat(' ',1,n);



fprintf(b,firstlines);
for i = 1:length(paramlines)
    fprintf(b,'  '); % Two spaces
    fprintf(b,paramlines{i}); % Dump the first column
    fprintf(b,repmat(' ',1,tabstop-length(paramlines{i})-2)); % pad spaces
    
    %Dump the second column, which has all sorts of formatting ugliness.
    if isstr(defaults{i})
        fprintf(b,repmat(' ',1,4)); % 4 spaces
        fprintf(b,defaults{i});    % This had better by 1 char (y/n) or things will get out of alignment. 
        if ~isempty(units{i})
            fprintf(b,repmat(' ',1,6)); % 6 spaces, but only if we need to add units.
        end
    elseif is_value_int(i)
        fprintf(b,repmat(' ',1,4)); % 4 spaces
        fprintf(b,num2str(defaults{i}));    % This had better by 1 char (y/n)
        if ~isempty(units{i})
            fprintf(b,repmat(' ',1,7-length(num2str(defaults{i})))); % 6 spaces, but written this way so a month of 11 or 12 won't produce too many spaces.
        end
    else
        fprintf(b,'%11.6f',defaults{i});
    end
    
    % Dump the third column, but only if it exists.
    
    if ~isempty(units{i})
        fprintf(b,repmat(' ',1,6));
        fprintf(b,units{i});
        
    end
    
    
    
    
    
    
    fprintf(b,'\n');
    
end

fclose(b);