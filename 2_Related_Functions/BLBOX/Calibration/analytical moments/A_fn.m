 function A = A_fn(kk,xx,eps,nu,alph,n)
      
      vec = 0:1:n;
      A = [];
      for j = 1:length(xx) 
           A(j) = 1- ((nu+xx(j))*eps)^(alph-kk)/gamma(alph-kk) * ...
                  sum((-1)^vec * ((nu+xx(j))*eps)^vec / ((alph-kk+vec)*factorial(vec)));
      end
     
   