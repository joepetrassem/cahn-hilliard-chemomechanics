function results = postProcessing(raw, params)

nbin = params.nbin;
rbin = params.r_bin';
redge = params.r_edge';
dr = params.dr;

c_ref = params.c_ref;
alph = params.alph;
nu = params.nu;
kappa = params.kappa;
gamm = params.gamm;
D0 = params.D0;
%R_gT = params.R_gT;
%Far = params.Far;

sol_total = raw.sol_total;
t_total = raw.t_total;
I = raw.I;

sol_u = sol_total(:,1:nbin);
sol_c = sol_total(:,nbin+1:2*nbin);
sol_ep = sol_total(:,2*nbin+1);


SOC = zeros(1,length(t_total));
for i=1:length(t_total)
    SOC(i) = trapz(rbin,3.*(sol_c(i,:).*rbin.^2));
end


%calculating sigma_kk, mu, F
sigma_kk = zeros(length(t_total),nbin);
sigma_rr = zeros(length(t_total),nbin);
sigma_tt = zeros(length(t_total),nbin);
e_elast_kk = zeros(length(t_total),nbin);
e_rr_el = zeros(length(t_total),nbin);
e_tt_el = zeros(length(t_total),nbin);
work_el = zeros(length(t_total),nbin);
mu = zeros(length(t_total),nbin);
F = zeros(length(t_total),nbin+1);
for i = 1:length(t_total)
    %get variables for time step
    w = sol_total(i,:);
    u = [0,w(1:nbin)];
    c = w(nbin+1:2*nbin);
    ept = w(end);

    %calculate strains
    e_s = c-c_ref; 

    e_rr = diff(u)./dr;

    e_tt = (u(2:end)+u(1:end-1))./(2.*rbin);

    e_rr_el(i,:) = e_rr - e_s;
    e_tt_el(i,:) = e_tt - e_s;

    e_elast_kk(i,:) = (e_rr + 2*e_tt - 3*e_s); %divide by 3 for volumetric

    %calculate stresses 
    sigma_rr(i,:) = alph*(e_rr - e_s + nu*(e_rr + 2*e_tt - 3*e_s));
    sigma_tt(i,:) = alph*(e_tt - e_s + nu*(e_rr + 2*e_tt - 3*e_s));

    %calculate volumetric stress
    sigma_kk(i,:) = (sigma_rr(i,:) + 2.*sigma_tt(i,:));   %divide by 3 for volumetric

    work_el(i,:) = sigma_rr(i,:).*e_rr_el(i,:) + 2*sigma_tt(i,:).*e_tt_el(i,:);

    %calulate mu
    J = [0,redge(2:end-1).^2.*diff(c)/dr, 0]; % dcdx with two boundary conditions enforced
    dJdx = diff(J)./(dr.*rbin.^2); % d2cdx2 at midpoints
    
    mu(i,:) = mu_chem(c) - kappa*dJdx - gamm.*sigma_kk(i,:)./3;
    %general function of c for now, need to infer what it looks like from data
    %and check scalings and parameters
    

    c_edge = (c(2:end) + c(1:end-1))/2; %interpolates for c on edges


    F(i,2:end-1) = -D0.*c_edge.*(1-c_edge).*(diff(mu(i,:))/dr);


    F(i,1) = 0; %Boundary conditions
    F(i,end) = I(i);
end

results.t_total = t_total;
results.SOC = SOC;
results.c = sol_c;
results.u = sol_u;
results.ep = sol_ep;
results.sigma_rr = sigma_rr;
results.sigma_tt = sigma_tt;
results.sigma_kk = sigma_kk;
results.I = I;
%results.e_elast_kk = e_elast_kk;
%results.e_rr_el = e_rr_el;
%results.e_tt_el = e_tt_el;
%results.work_el = work_el;
results.mu = mu;
%results.F =  F;


end