% # useful function for expectations of eta for truncated gamma
% # n indicates how many terms we include in the summation which approximates the integral 
% # epsilon is the value which we integrate from i.e. the Gamma distribution for eta has support
% # (epsilon, inf)

function T=Exp_fn_e(k,x,eps, nu, alpha, n) 
if nargin<6
n=100;
end
alph=alpha;
T=zeros(size(alph));
ind=find(alpha==1);
alph(ind)= 1.000001;
ind=find(alpha==2);
alph(ind)=2.000001;
ind=find(alpha==3);
alph(ind)=3.000001;
ind=find(alpha==4);
alph(ind)=4.000001;
     
   

   if (eps == 0)
       for i=1:length(alph)           
           if numel(x)==1
               T(i)=(nu(i)/(nu(i)+x))^alph(i) * (nu(i)+x)^k * prod(alph(i)-1:-1:alph(i)-k)^(-1); 
           else
               T(i)=(nu(i)/(nu(i)+x(i)))^alph(i) * (nu(i)+x(i))^k * prod(alph(i)-1:-1:alph(i)-k)^(-1); 
           end
       end
   
   else   
       sum1=zeros(length(alph),1);
       sum2=sum1;
       
    
       j=0:n;    
       N = j(:);
       j(N>170) = 171;
       m = max([1; j(:)]);
       N = [1 1 cumprod(2:m)];
       fac = N(j+1);
       
       
           sum1=sum((-1).^j.*(x+nu).^j.*eps.^j./((alph-k+j).*fac),2);
           sum2=sum((-1).^j.*nu.^j.*eps.^j./((alph+j).*fac),2);
           T=(gamma(alph-k)./(x+nu).^(alph-k) - ...
               eps.^(alph-k).*sum1)./(gamma(alph)./nu.^(alph) - ...
               eps.^alph.*sum2);
    end
           
   end
  
