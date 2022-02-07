function label=seasonname(nr,language)

if strcmp(language,'EN')
    seasons={'Winter','Spring','Summer','Fall'};
elseif strcmp(language,'NL')
    seasons={'Winter','Lente','Zomer','Herfst'};
end

label=seasons{nr};