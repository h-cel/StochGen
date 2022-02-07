function z = h(family,x,y,d,condvar,type)
%% h Implementation of the h-function of Aas et al. (2009)
%   function for calculations with (inverse) conditional copula C(u,v).
%   This function was first proposed by Aas et al. (2009) and further
%   expanded in subsequent papers and books.
%
%   Inputs:
%       family: copula family (string). Currently implemented: Frank,
%       Clayton, Gumbel, Gaussian and t
%       type = 'norm'
%           x and y: u and v values on [0,1]
%       type = 'inv'
%           x: probability level of conditional copula
%           y: value of u or v, depending on conditioning variable
%       d: copula parameter
%       condvar: conditioning variable: 'U' or 'V'
%       type: 'norm' or 'inv', type of h-function needed
%   Output:
%       z: a probability level in case of 'norm'
%          a value for u or v in case of 'inv'
%
%   Originally implemented by M. T. Pham (?) and H. Vernieuwe
%   Last update by J. Van de Velde on 14/07/'21: fixed bugs

%% Error check

if ~isreal(x) && imag(x) < 0.0001
    x = real(x);
elseif ~isreal(x) && imag(x) > 0.0001
    z = NaN;
    return
end

if ~isreal(y) && imag(y) < 0.0001
    y = real(x);
elseif ~isreal(y) && imag(y) > 0.0001
    z = NaN;
    return
end

%% Implementation

if d==0
    d=0.0001;
end

switch family
    case 'Frank'
        switch type
            case 'norm'
                u=x;
                v=y;
                switch condvar
                    case 'U'
                        z = -(exp(-d.*u).*(exp(-d.*v)-1)./(-exp(-d)...
                            -exp(-d.*u-d.*v)+exp(-d.*u)+exp(-d.*v)));
                    case 'V'
                        z = -(exp(-(d.*u))-1).*exp(-(d.*v))./(-exp(-d)...
                            -exp((-d.*u-d.*v))+exp(-(d.*u))+exp(-(d.* v)));
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    
                    prec=1000;
                    iter=0;
                    
                    left=0;
                    right=1;
                    
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=y(i);
                                v2=midpoint;
                                k = -(exp(-d.*u).*(exp(-d.*v2)-1)./(-exp(-d)-...
                                    exp(-d.*u-d.*v2)+exp(-d.*u)+exp(-d.*v2)));
                            case 'V'
                                v=y(i);
                                u2=midpoint;
                                k = -(exp(-(d.*u2))-1).*exp(-(d.*v))...
                                    ./(-exp(-d)-exp((-d.*u2-d.*v))...
                                    +exp(-(d.*u2))+exp(-(d.* v)));
                        end
                        if k > t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end
                    
                    z(i)=(right+left)/2;
                end
        end
    case 'Clayton' % Implementation based on Aas et al. (2009)
        switch type
            case 'norm'
                u=x;
                v=y;
                switch condvar
                    case 'U'
                        z= u.^(-d-1).*(v.^(-d)+u.^(-d)-1).^(-1-1/d);
                    case 'V'
                        z= v.^(-d-1).*(u.^(-d)+v.^(-d)-1).^(-1-1/d);
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    prec=1000;
                    iter=0;
                    left=0;
                    right=1;
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=y(i);
                                v2=midpoint;
                                k = u.^(-d-1).*(v2.^(-d)+u.^(-d)-1).^(-1-1/d);
                            case 'V'
                                v=y(i);
                                u2=midpoint;
                                k = v.^(-d-1).*(u2.^(-d)+v.^(-d)-1).^(-1-1/d);
                        end
                        if k > t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end
                    
                    z(i)=(right+left)/2;
                end
        end
    case 'Clayton 90°'
        switch type
            case 'norm'
                switch condvar
                    case 'U'
                        u=x;
                        v=1-y;
                        z= 1- (u.^(-d-1).*(v.^(-d)+u.^(-d)-1).^(-1-1/d));
                    case 'V'
                        u=1-x;
                        v=y;
                        z= 1-( v.^(-d-1).*(u.^(-d)+v.^(-d)-1).^(-1-1/d));
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    prec=1000;
                    iter=0;
                    left=0;
                    right=1;
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=y(i);
                                v2=1-midpoint;
                                k =1-(u.^(-d-1).*(v2.^(-d)+u.^(-d)-1).^(-1-1/d));
                                
                            case 'V'
                                v=y(i);
                                u2=1-midpoint;
                                k = 1-(v.^(-d-1).*(u2.^(-d)+v.^(-d)-1).^(-1-1/d));
                                
                        end
                        if k > t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end
                    
                    z(i)=(right+left)/2;

                end
        end
    case 'Clayton 180°' %Survival Clayton
        switch type
            case 'norm'
                u=1-x;
                v=1-y;
                switch condvar
                    case 'U'
                        z= 1- (u.^(-d-1).*(v.^(-d)+u.^(-d)-1).^(-1-1/d));
                    case 'V'
                        z= 1- (v.^(-d-1).*(u.^(-d)+v.^(-d)-1).^(-1-1/d));
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    prec=1000;
                    iter=0;
                    left=0;
                    right=1;
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=1-y(i);
                                v2=1-midpoint;
                                k =1- (u.^(-d-1).*(v2.^(-d)+u.^(-d)-1).^(-1-1/d));
                            case 'V'
                                v=1-y(i);
                                u2=1-midpoint;
                                k = 1- (v.^(-d-1).*(u2.^(-d)+v.^(-d)-1).^(-1-1/d));
                        end
                        if k > t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end
                    
                    z(i)=(right+left)/2;
                end
        end
    case 'Clayton 270°'
        switch type
            case 'norm'
                switch condvar
                    case 'U'
                        u=1-x;
                        v=y;
                        z= u.^(-d-1).*(v.^(-d)+u.^(-d)-1).^(-1-1/d);
                    case 'V'
                        u=x;
                        v=1-y;
                        z= v.^(-d-1).*(u.^(-d)+v.^(-d)-1).^(-1-1/d);
                end
                if ~isreal(z) && imag(z) < 0.00001 %Arbitrary threshold
                    z = real(z);
                elseif real((v.^(-d)+u.^(-d)-1).^(-1-1/d)) < 0.00001
                    z = u.^(-d-1).*0;
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    prec=1000;
                    iter=0;
                    left=0;
                    right=1;
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=1-y(i);
                                v2=midpoint;
                                k = u.^(-d-1).*(v2.^(-d)+u.^(-d)-1).^(-1-1/d);
                            case 'V'
                                v=1-y(i);
                                u2=midpoint;
                                k = v.^(-d-1).*(u2.^(-d)+v.^(-d)-1).^(-1-1/d);
                                
                        end
                        if k > t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end

                    
                    z(i)=(right+left)/2;
                end
        end
    case 'Gumbel'
        switch type
            case 'norm'
                u=x;
                v=y;
                switch condvar
                    case 'U'
                        z = -((-log(u)).^d+(-log(v)).^d).^(1./d).*(-log(u)).^d./u./log(u)./...
                            ((-log(u)).^d+(-log(v)).^d).*exp(-((-log(u)).^d+(-log(v)).^d).^(1./d));
                    case 'V'
                        z = -((-log(u)).^d+(-log(v)).^d).^(1./d).*(-log(v)).^d./v./log(v)./...
                            ((-log(u)).^d+(-log(v)).^d).*exp(-((-log(u)).^d+(-log(v)).^d).^(1./d));
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    prec=1000;
                    iter=0;
                    left=0;
                    right=1;
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=y(i);
                                v2=midpoint;
                                k = -((-log(u)).^d+(-log(v2)).^d)^(1./d)*(-log(u)).^d./u./log(u)./...
                                    ((-log(u)).^d+(-log(v2)).^d).*exp(-((-log(u)).^d+(-log(v2)).^d)^(1./d));
                            case 'V'
                                v=y(i);
                                u2=midpoint;
                                k = -((-log(u2)).^d+(-log(v)).^d).^(1./d).*(-log(v)).^d./v./log(v)./...
                                    ((-log(u2)).^d+(-log(v)).^d).*exp(-((-log(u2)).^d+(-log(v)).^d)^(1./d));
                        end
                        if k > t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end
                    
                    z(i)=(right+left)/2;
                end
        end
    case 'Gumbel 90°'
        switch type
            case 'norm'
                switch condvar
                    case 'U'
                        u=x;
                        v=1-y;
                        z =1-( -((-log(u)).^d+(-log(v)).^d).^(1./d).*(-log(u)).^d./u./log(u)./...
                            ((-log(u)).^d+(-log(v)).^d).*exp(-((-log(u)).^d+(-log(v)).^d).^(1./d)));
                    case 'V'
                        u=1-x;
                        v=y;
                        z =1-( -((-log(u)).^d+(-log(v)).^d).^(1./d).*(-log(v)).^d./v./log(v)./...
                            ((-log(u)).^d+(-log(v)).^d).*exp(-((-log(u)).^d+(-log(v)).^d).^(1./d)));
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    prec=1000;
                    iter=0;
                    left=0;
                    right=1;
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=y(i);
                                v2=1-midpoint;
                                k =1-( -((-log(u)).^d+(-log(v2)).^d)^(1./d)*(-log(u)).^d./u./log(u)./...
                                    ((-log(u)).^d+(-log(v2)).^d).*exp(-((-log(u)).^d+(-log(v2)).^d)^(1./d)));
                                
                                
                            case 'V'
                                v=y(i);
                                u2=1-midpoint;
                                k = 1-(-((-log(u2)).^d+(-log(v)).^d).^(1./d).*(-log(v)).^d./v./log(v)./...
                                    ((-log(u2)).^d+(-log(v)).^d).*exp(-((-log(u2)).^d+(-log(v)).^d)^(1./d)));
                                
                        end
                        if k > t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end
                    z(i)=(right+left)/2;
                end
        end
    case 'Gumbel 180°' %survival Gumbel
        switch type
            case 'norm'
                u=1-x;
                v=1-y;
                switch condvar
                    case 'U'
                        z =1-( -((-log(u)).^d+(-log(v)).^d).^(1./d).*(-log(u)).^d./u./log(u)./...
                            ((-log(u)).^d+(-log(v)).^d).*exp(-((-log(u)).^d+(-log(v)).^d).^(1./d)));
                    case 'V'
                        z =1-( -((-log(u)).^d+(-log(v)).^d).^(1./d).*(-log(v)).^d./v./log(v)./...
                            ((-log(u)).^d+(-log(v)).^d).*exp(-((-log(u)).^d+(-log(v)).^d).^(1./d)));
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    prec=1000;
                    iter=0;
                    left=0;
                    right=1;
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=1-y(i);
                                v2=1-midpoint;
                                k =1-( -((-log(u)).^d+(-log(v2)).^d)^(1./d)*(-log(u)).^d./u./log(u)./...
                                    ((-log(u)).^d+(-log(v2)).^d).*exp(-((-log(u)).^d+(-log(v2)).^d)^(1./d)));
                            case 'V'
                                v=1-y(i);
                                u2=1-midpoint;
                                k = 1-(-((-log(u2)).^d+(-log(v)).^d).^(1./d).*(-log(v)).^d./v./log(v)./...
                                    ((-log(u2)).^d+(-log(v)).^d).*exp(-((-log(u2)).^d+(-log(v)).^d)^(1./d)));
                        end
                        if k > t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end
                    
                    z(i)=(right+left)/2;
                end
        end
    case 'Gumbel 270°'
        switch type
            case 'norm'
                switch condvar
                    case 'U'
                        u=1-x;
                        v=y;
                        z = -((-log(u)).^d+(-log(v)).^d).^(1./d).*(-log(u)).^d./u./log(u)./...
                            ((-log(u)).^d+(-log(v)).^d).*exp(-((-log(u)).^d+(-log(v)).^d).^(1./d));
                    case 'V'
                        u=x;
                        v=1-y;
                        z = -((-log(u)).^d+(-log(v)).^d).^(1./d).*(-log(v)).^d./v./log(v)./...
                            ((-log(u)).^d+(-log(v)).^d).*exp(-((-log(u)).^d+(-log(v)).^d).^(1./d));
                end
            case 'inv'
                t2=x;
                z=zeros(size(x));
                for i=1:length(x)
                    prec=1000;
                    iter=0;
                    left=0;
                    right=1;
                    while prec > 0.0000001
                        midpoint=(right+left)/2;
                        switch condvar
                            case 'U'
                                u=1-y(i);
                                v2=midpoint;
                                k = -((-log(u)).^d+(-log(v2)).^d)^(1./d)*(-log(u)).^d./u./log(u)./...
                                    ((-log(u)).^d+(-log(v2)).^d).*exp(-((-log(u)).^d+(-log(v2)).^d)^(1./d));
                            case 'V'
                                v=1-y(i);
                                u2=midpoint;
                                k = -((-log(u2)).^d+(-log(v)).^d).^(1./d).*(-log(v)).^d./v./log(v)./...
                                    ((-log(u2)).^d+(-log(v)).^d).*exp(-((-log(u2)).^d+(-log(v)).^d)^(1./d));
                                
                        end
                        if k < t2(i)
                            right=midpoint;
                        else
                            left=midpoint;
                        end
                        iter=iter+1;
                        prec=(1/2)^iter;
                    end
                    z(i)=(right+left)/2;
                end
        end
    case 'Gaussian'
        % formules voor conditionele in Meyer C*, 2013,
        % bivariate standard normal, conditional Y|X is again normally
        % distributed with mean rho*x and standard deviation sqrt(1-rho^2)
        % numerically checked ((C(u,v+deltav)-C(u,v))/deltav)
        switch type
            case 'norm'
                u=x;
                v=y;
                switch condvar
                    case 'U'
                        z=normcdf((norminv(v)-d*norminv(u))/(sqrt(1-d^2)));
                    case 'V'
                        z=normcdf((norminv(u)-d*norminv(v))/(sqrt(1-d^2)));
                end
            case 'inv'
                t2=x;
                switch condvar
                    case 'U'
                        u=y;
                        z=normcdf(norminv(t2).*sqrt(1-d^2)+d*norminv(u));
                    case 'V'
                        v=y;
                        z=normcdf(norminv(t2).*sqrt(1-d^2)+d*norminv(v));
                end
        end
    case 't'
        switch type
            case 'norm'
                u=x;
                v=y;
                switch condvar
                    case 'U'
                        z=tcdf((tinv(v,d(2))-d(1)*tinv(u,d(2)))./sqrt((d(2)+tinv(u,d(2)).^2)*(1-d(1)^2)/(d(2)+1)),d(2)+1);
                    case 'V'
                        z=tcdf((tinv(u,d(2))-d(1)*tinv(v,d(2)))./sqrt((d(2)+tinv(v,d(2)).^2)*(1-d(1)^2)/(d(2)+1)),d(2)+1);
                end
            case 'inv'
                t2=x;
                switch condvar
                    case 'U'
                        u=y;
                        z=tcdf(tinv(t2,d(2)+1).*sqrt(((d(2)+tinv(u,d(2)).^2).*(1-d(1)^2))/(d(2)+1))+d(1).*tinv(u,d(2)),d(2));
                        
                    case 'V'
                        v=y;
                        z=tcdf(tinv(t2,d(2)+1).*sqrt(((d(2)+tinv(v,d(2)).^2).*(1-d(1)^2))/(d(2)+1))+d(1).*tinv(v,d(2)),d(2));
                end
        end
end
end
