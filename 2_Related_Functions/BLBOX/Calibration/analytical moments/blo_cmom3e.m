function thirdmoment=blo_cmom3e(par,h)

lambda=par(:,1);
beta=par(:,2);
gamma=par(:,3);
mux=par(:,4);
eta=par(:,5);

phi=gamma./eta;
kappa=beta./eta;
muc=1+kappa./phi;

mux2=2*mux.^2;
mux3=6*mux.^3;



mean=h*lambda/eta*muc*mux;
var=2*lambda/eta*muc*(mux2+beta*mux^2/gamma)*h/eta ...
    +2*lambda/eta*muc*mux^2*beta*eta*(1-exp(-gamma*h))...
    /(gamma^2*(gamma^2-eta^2)) ...
    -2*lambda/eta*muc*(mux2+beta*gamma*mux^2/(gamma^2-eta^2))...
    *(1-exp(-eta*h))/(eta^2);

T=6*((1/6.*(12.*(phi*eta)^9.*mux3*exp(h*(eta))*exp(h*(phi*eta))+ ...
    6.*(eta)^4 ...
    *mux2*mux*(phi*eta)^4.*((kappa*eta))*exp(h*(eta))^2.+24.*(eta)^5.* ...
    (phi*eta)^4.* ...
    mux2*mux*h*((kappa*eta))*exp(h*(eta))*exp(h*(phi*eta))-36.*(eta)^3.* ...
    mux2* ...
    mux*(phi*eta)^7*h^2*(lambda)*(muc)*exp(h*(eta))^2*exp(h*(phi*eta))-30.* ...
    (eta)^3*mux2*mux*(phi*eta)^6*h*((kappa*eta))*exp(h*(eta))*exp(h*(phi*eta))+12 ...
    *(eta)^4*(phi*eta)^3*mux^3*((kappa*eta))^2*exp(h*(eta))^2+12.*(phi*eta)^7* ...
    mux^3*((kappa*eta))^2*exp(h*(eta))*exp(h*(phi*eta))+6.*(eta)*(phi*eta)^9*h* ...
    mux3*exp(h*(eta))*exp(h*(phi*eta))-48.*(eta)^7*mux^3*((kappa*eta))^2* ...
    exp(h*(eta))^2-12.*(phi*eta)^9*exp(h*(eta))^2*exp(h*(phi*eta))*mux3 ...
    +6.*(phi*eta)^8*mux^3*h*((kappa*eta))*(lambda)*(muc)*exp(h*(eta))*exp(h* ...
    (phi*eta))+(eta)*mux^3*(phi*eta)^9*(lambda)^2*(muc)^2*h^3*exp(h*(eta ...
    ))^2*exp(h*(phi*eta))-4.*(eta)^7*(phi*eta)^3*mux^3*(lambda)^2* ...
    (muc)^2 ...
    *h^3*exp(h*(eta))^2*exp(h*(phi*eta))-30.*(eta)^2*(phi*eta)^6*mux^3*h* ...
    ((kappa*eta))*(lambda)*(muc) ...
    *exp(h*(eta))*exp(h*(phi*eta))-24.*(eta)^6*mux2 ...
    *mux*((kappa*eta))*(phi*eta)^2*exp(h*(eta))*exp(h*(phi*eta))+6.*(eta)^2* ...
    (phi*eta)^5 ...
    *mux^3*((kappa*eta))^2*exp(h*(eta))^2-24.*(eta)^7*(phi*eta)^3*mux2*mux ...
    *h^2*(lambda)*(muc)*exp(h*(eta))^2*exp(h*(phi*eta))-36.*(eta)^3* ...
    (phi*eta)^6*mux^3*h^2*((kappa*eta))*(lambda)*(muc)*exp(h*(eta))^2*exp(h* ...
    (phi*eta))-24.*(eta)^7*mux^3*h^2*((kappa*eta))*(lambda)*(muc)*(phi*eta)^2* ...
    exp(h*(eta))^2*exp(h*(phi*eta))+6.*(eta)*(phi*eta)^8*mux^3*h^2*((kappa*eta)) ...
    *(lambda)*(muc)*exp(h*(eta))^2*exp(h*(phi*eta))-18.*(eta)^3*(phi*eta)^4 ...
    *mux^3*((kappa*eta))^2*exp(h*(eta))*exp(h*(phi*eta))-72.*(eta)^2*(phi*eta)^7* ...
    mux3*exp(h*(eta))*exp(h*(phi*eta))-6.*(eta)^3*(phi*eta)^5*mux^3*h*( ...
    (kappa*eta))^2*exp(h*(eta))^2+84.*(eta)^5*mux^3*((kappa*eta))^2*(phi*eta)^2* ...
    exp(h*(eta))^2-6.*(eta)^3*mux2*mux*(phi*eta)^5*((kappa*eta))*exp(h*(eta) ...
    )*exp(h*(phi*eta))-48.*(eta)^6*mux3*(phi*eta)^3*exp(h*(eta))*exp(h* ...
    (phi*eta))+30.*(eta)^5* ...
    mux^3*h*((kappa*eta))*(lambda)*(muc)*(phi*eta)^3*exp(h*( ...
    eta))^2-36.*(eta)^2*mux2*mux*(phi*eta)^7*h*(lambda)*(muc)*exp(h*( ...
    eta))*exp(h*(phi*eta))-48.*(eta)^7*mux2*mux*((kappa*eta))*(phi*eta)*exp(h*(eta ...
    ))^2+6.*(eta)*mux2*mux*(phi*eta)^8*h*((kappa*eta))*exp(h*(eta))*exp(h* ...
    (phi*eta))-24.*(eta)^7*mux^3*h*((kappa*eta))^2*(phi*eta)*exp(h*(eta))^2+30.*( ...
    eta)^5*mux^3*h*((kappa*eta))^2*(phi*eta)^3*exp(h*(eta))^2+54.*(eta)^5* ...
    (phi*eta)^4*mux^3*h^2*((kappa*eta))*(lambda)*(muc)*exp(h*(eta))^2*exp(h* ...
    (phi*eta))-24.*(eta)^6*mux2*mux*((kappa*eta))*(phi*eta)^2*exp(h*(eta))^2+ ...
    48.* ...
    exp(h*(eta))^2*exp(h*(phi*eta))*mux^3*((kappa*eta))^2*(eta)^7+12. ...
    *(phi*eta)^8 ...
    *exp(h*(eta))^2*exp(h*(phi*eta))*h*mux*mux2*((kappa*eta))*(eta)-6. ...
    *(phi*eta)^9 ...
    *exp(h*(eta))^2*exp(h*(phi*eta))*mux2*mux*h*(lambda)*(muc)+6.* ...
    (phi*eta)^9*exp(h*(eta))^2*exp(h*(phi*eta))*h*mux3*(eta)-21. ...
    *(phi*eta)^8 ...
    *exp(h*(eta))^2*exp(h*(phi*eta))*mux2*mux*((kappa*eta))-24. ...
    *(phi*eta)^4*exp(h* ...
    (eta))^2*exp(h*(phi*eta))*mux^3*h*((kappa*eta))*(lambda)*(muc)*(eta)^4- ...
    138. ...
    *(phi*eta)^4*exp(h*(eta))^2*exp(h*(phi*eta))*mux*mux2*((kappa*eta))*(eta)^4+ ...
    108.*(phi*eta)^4*exp(h*(eta))^2*exp(h*(phi*eta))*mux*mux2*h*((kappa*eta))*( ...
    eta)^5+18.*(phi*eta)^4*exp(h*(eta))^2*exp(h*(phi*eta))*mux^3*((kappa*eta))^2 ...
    *(eta)^3+39.*(phi*eta)^5*exp(h*(eta))^2*exp(h*(phi*eta))*mux^3*( ...
    (kappa*eta))^2*(eta)^2+54.*(eta)^5*mux2*mux*(phi*eta)^5*h^2*(lambda) ...
    *(muc)*exp(h*(eta))^2*exp(h*(phi*eta))-36.*(eta)^3*(phi*eta)^7*h* ...
    mux3*exp(h*(eta))*exp(h*(phi*eta))+6.*(phi*eta)^9*h*mux*mux2*( ...
    lambda)*(muc)*exp(h*(eta))*exp(h*(phi*eta))+24.*(eta)^5*mux2*mux* ...
    ((kappa*eta))*(phi*eta)^3*exp(h*(eta))*exp(h*(phi*eta))-6.*(eta)^3*(phi*eta)^7* ...
    mux^3*(lambda)^2*(muc)^2*h^3*exp(h*(eta))^2*exp(h*(phi*eta))- ...
    18.*(eta)^3*(phi*eta)^4*mux^3*((kappa*eta))^2*exp(h*(eta))^2-24. ...
    *(eta)^7 ...
    *h*mux3*(phi*eta)^3*exp(h*(eta))*exp(h*(phi*eta))+6.*(eta)*mux2*mux ...
    *(phi*eta)^9*h^2*(lambda)*(muc)*exp(h*(eta))^2*exp(h*(phi*eta))+108.* ...
    (eta)^4*(phi*eta)^5*mux3*exp(h*(eta))*exp(h*(phi*eta))-132.*(eta)^2 ...
    *mux2*mux*(phi*eta)^6*((kappa*eta))*exp(h*(eta))*exp(h*(phi*eta))+24.*(eta)^4 ...
    *mux^3*h*((kappa*eta))*(lambda)*(muc)*(phi*eta)^4*exp(h*(eta))*exp(h*(phi*eta)) ...
    +9*(eta)^5*(phi*eta)^5*mux^3*(lambda)^2*(muc)^2*h^3*exp(h*( ...
    eta))^2*exp(h*(phi*eta))+36.*(eta)^5*mux2*mux*((kappa*eta))*(phi*eta)^3* ...
    exp(h*(eta))^2+54.*(eta)^5*(phi*eta)^5*h*mux3*exp(h*(eta))* ...
    exp(h*(phi*eta))-24.*(eta)^6*(phi*eta)^3*mux2*mux*h*(lambda)*(muc)* ...
    exp(h*(eta))*exp(h*(phi*eta))+150.*(eta)^4 ...
    *mux2*mux*(phi*eta)^4*((kappa*eta))* ...
    exp(h*(eta))*exp(h*(phi*eta))+12.*(eta)^4*(phi*eta)^3*mux^3*((kappa*eta))^2* ...
    exp(h*(eta))*exp(h*(phi*eta))-6.*(eta)^3*mux2*mux*(phi*eta)^5*((kappa*eta))* ...
    exp(h*(eta))^2-6.*(eta)^3*(phi*eta)^5*mux^3*h*((kappa*eta))*(lambda)*(muc ...
    )*exp(h*(eta))^2+24.*mux2*mux*(phi*eta)^8*((kappa*eta))*exp(h*(eta)) ...
    *exp(h*(phi*eta))-24.*(eta)^7*mux^3*h*((kappa*eta))*(lambda)*(muc)*(phi*eta)* ...
    exp(h*(eta))^2+54.*(eta)^4*mux2*mux*(phi*eta)^5*h*(lambda)*(muc) ...
    *exp(h*(eta))*exp(h*(phi*eta))-42.*(eta)^2*(phi*eta)^5*mux^3*((kappa*eta))^2 ...
    *exp(h*(eta))*exp(h*(phi*eta))+24.*(phi*eta)*exp(h*(eta))^2*exp(h*(phi*eta))* ...
    mux^3*h*((kappa*eta))*(lambda)*(eta)^7*(muc)-48.*(phi*eta)^2* ...
    exp(h*(eta))^2*exp(h*(phi*eta))*mux*mux2*h*((kappa*eta))*(eta)^7+48.* ...
    (phi*eta)*exp(h*(eta))^2*exp(h*(phi*eta))*mux*mux2*((kappa*eta))*(eta)^7-24. ...
    *(phi*eta)*exp(h*(eta))^2*exp(h*(phi*eta))*mux^3*h*((kappa*eta))^2*(eta)^7+ ...
    24.*(phi*eta)^2*exp(h*(eta))^2*exp(h*(phi*eta))*mux*mux2*((kappa*eta))*( ...
    eta)^6+72.*(phi*eta)^7*exp(h*(eta))^2*exp(h*(phi*eta))*mux3* ...
    (eta)^2-6.*(phi*eta)^8*exp(h*(eta))^2*exp(h*(phi*eta))*mux^3*h*((kappa*eta) ...
    )*(lambda)*(muc)-30.*(phi*eta)^3*exp(h*(eta))^2*exp(h*(phi*eta))*mux^3 ...
    *h*((kappa*eta))*(lambda)*(eta)^5*(muc)-36.*(phi*eta)^3*exp(h*(eta))^2* ...
    exp(h*(phi*eta))*mux*mux2*((kappa*eta))*(eta)^5-12.*(phi*eta)^3* ...
    exp(h*(eta))^2*exp(h*(phi*eta))*mux^3*((kappa*eta))^2*(eta)^4+54.* ...
    (phi*eta)^3*exp(h*(eta))^2*exp(h*(phi*eta))*mux^3*h*((kappa*eta))^2*(eta)^5 ...
    -24.*(phi*eta)^3*exp(h*(eta))^2*exp(h*(phi*eta))*h*mux3*(eta)^7+24. ...
    *(phi*eta)^3*exp(h*(eta))^2*exp(h*(phi*eta))*mux2*mux*h*(lambda)*(muc ...
    )*(eta)^6-84.*(phi*eta)^2*exp(h*(eta))^2*exp(h*(phi*eta))*mux^3*( ...
    (kappa*eta))^2*(eta)^5+54.*(phi*eta)^5*exp(h*(eta))^2*exp(h*(phi*eta))*h* ...
    mux3*(eta)^5-54.*(phi*eta)^5*exp(h*(eta))^2*exp(h*(phi*eta))* ...
    mux2*mux*h*(lambda)*(muc)*(eta)^4+6.*(phi*eta)^5*exp(h*(eta))^2 ...
    *exp(h*(phi*eta))*mux^3*h*((kappa*eta))*(lambda)*(eta)^3*(muc)-108.* ...
    (phi*eta)^5*exp(h*(eta))^2*exp(h*(phi*eta))*mux3*(eta)^4+6.* ...
    (phi*eta)^5 ...
    *exp(h*(eta))^2*exp(h*(phi*eta))*mux*mux2*((kappa*eta))*(eta)^3-72.* ...
    (phi*eta)^6*exp(h*(eta))^2*exp(h*(phi*eta))*mux*mux2*h*((kappa*eta))* ...
    (eta)^3+117.*(phi*eta)^6*exp(h*(eta))^2*exp(h*(phi*eta))*mux*mux2 ...
    *((kappa*eta))*(eta)^2+30.*(phi*eta)^6*exp(h*(eta))^2*exp(h*(phi*eta))*mux^3 ...
    *h*((kappa*eta))*(lambda)*(muc)*(eta)^2+6.*(phi*eta)^7*exp(h*(eta))^2* ...
    exp(h*(phi*eta))*h*mux^3*((kappa*eta))^2*(eta)-9.*(phi*eta)^7*exp(h*(eta))^2 ...
    *exp(h*(phi*eta))*mux^3*((kappa*eta))^2-36.*(phi*eta)^7*exp(h*(eta))^2* ...
    exp(h*(phi*eta))*h*mux3*(eta)^3-36.*(phi*eta)^5*exp(h*(eta))^2  ...
    *exp(h*(phi*eta))*mux^3*h*((kappa*eta))^2*(eta)^3+36.*(phi*eta)^7* ...
    exp(h*(eta))^2*exp(h*(phi*eta))*mux2*mux*h*(lambda)*(muc)* ...
    (eta)^2+48.*(phi*eta)^3*exp(h*(eta))^2*exp(h*(phi*eta))*mux3* ...
    (eta)^6-24.*(eta)^5*mux2*mux*((kappa*eta))*(phi*eta)^3*exp(h*(eta)) ...
    -12.*(eta)^4*(phi*eta)^3*mux^3*((kappa*eta))^2*exp(h*(eta))-12*(eta)^4 ...
    *mux2*mux*(phi*eta)^4*((kappa*eta))*exp(h*(phi*eta))-6.*(eta)^4*mux2*mux ...
    *(phi*eta)^4*((kappa*eta))*exp(h*(eta))+18.*(eta)^3*(phi*eta)^4*mux^3*( ...
    (kappa*eta))^2*exp(h*(eta))+6.*(eta)^3*(phi*eta)^5*mux2*mux*((kappa*eta))* ...
    exp(h*(eta))-3.*mux2*mux*(phi*eta)^8*((kappa*eta))*exp(h*(phi*eta))+24.*( ...
    eta)^6*mux2*mux*((kappa*eta))*(phi*eta)^2*exp(h*(eta))-3.*(phi*eta)^7* ...
    mux^3*((kappa*eta))^2*exp(h*(phi*eta))-6.*(eta)^2*(phi*eta)^5*mux^3*( ...
    (kappa*eta))^2*exp(h*(eta))+3.*(eta)^2*(phi*eta)^5*mux^3*((kappa*eta))^2 ...
    *exp(h*(phi*eta))+15.*(eta)^2*mux2*mux*(phi*eta)^6*((kappa*eta))*exp(h*(phi*eta) ...
    ))*lambda*muc/((eta)^2+2.*(phi*eta)*(eta)+(phi*eta)^2)/((phi*eta)^4-2.*( ...
    eta)*(phi*eta)^3-3.*(eta)^2*(phi*eta)^2+8.*(eta)^3*(phi*eta)-4.*(eta)^4) ...
    /exp(h*(eta))^2/(eta)^4/exp(h*(phi*eta))/(phi*eta)^3));

thirdmoment=T-3*mean*(var+mean^2)+2*mean^3;
