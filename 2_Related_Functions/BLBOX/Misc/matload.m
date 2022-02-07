%functie voor het laden van een mat-file.

function X=matload(matfile)

A=load(matfile);
vars=whos('-file',matfile);
X=A.(vars.name);

