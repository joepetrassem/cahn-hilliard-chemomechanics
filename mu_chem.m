function mu = mu_chem(c)
%MU_CHEM chemical component of potential mu
mu = log(c./(1-c)) + 3.*(1-2.*c);% - 0.5*c.^3;   %- exp(-((c-0.75)/0.1).^2); % + c.^3 for asymmetric wells.
%mu = log(c./(1-c)) - 2*exp(-0.5*((c - 0.8)/0.07).^2);
%mu = log(c./(1-c))  - 3*exp(-((c-0.91)./0.07).^2);

end

