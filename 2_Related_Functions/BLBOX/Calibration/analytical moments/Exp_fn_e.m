% # useful function for expectations of eta for truncated gamma
% # n indicates how many terms we include in the summation which approximates the integral 
% # epsilon is the value which we integrate from i.e. the Gamma distribution for eta has support
% # (epsilon, inf)

function T=Exp_fn_e(k,a,eps, nu, alpha) 
     
A = gammainc(alpha-k , (nu+a)*eps)*gamma(alpha-k)/(nu+a)^(alpha-k);

B = gammainc(alpha , nu* eps) *gamma(alpha) / nu^alpha;

T=A/B;

           
   end
  
