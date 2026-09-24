function [value, isterminal, direction] = elastEvent(t,w,wdot,params)
nbin = params.nbin;
u_surf = w(nbin);
c_surf = w(2*nbin);
alph = params.alph;
ept = w(end);
sig_y = params.sig_y;
A = params.A;

value = [abs(A*(alph*u_surf - ept))-sig_y]; %min(c_surf-0.01,0.99-c_surf)];

isterminal = [1]; %1];

direction = [1]; %0];
end