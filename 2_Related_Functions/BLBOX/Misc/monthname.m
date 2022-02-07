function label=monthname(nr,language)

if strcmp(language,'EN')
    months={'January','February','March','April','May','June','July', ...
        'August','September','October','November','December'};
elseif strcmp(language,'NL')
    months={'Januari','Februari','Maart','April','Mei','Juni','Juli', ...
        'Augustus','September','Oktober','November','December'};
end

label=months{nr};
