%% Function to create a realization of a poisson process with rate r

function [x]=randexp(r)

U=rand(1);
x=-log(U)./r;

end


