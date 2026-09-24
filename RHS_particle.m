function f = RHS_particle(t, w, wdot, params, I, regime)

nbin = params.nbin;
f = zeros(length(w),1);

%state vector
%w = [u(1:nbin) , c(1:nbin), ep]
u_in = w(1:nbin);
c = w(nbin+1:2*nbin);
ep=w(2*nbin+1);

c_dot = wdot(nbin+1:2*nbin);
ep_dot = wdot(2*nbin+1);


u = [0; u_in]; %boundary condition, fixed u at r=0


%% Particle mechanics

e_s = c-params.c_ref; 

e_rr = diff(u)./params.dr;

e_tt = (u(2:end)+u(1:end-1))./(2*params.r_bin);

sigma_rr = zeros(nbin+1,1); %ghost bin

sigma_rr(1:nbin) = params.alph.*(e_rr - e_s + params.nu*(e_rr + 2*e_tt - 3*e_s));  

sigma_tt = params.alph*(e_tt-e_s + params.nu*(e_rr + 2*e_tt - 3*e_s));

%% Surface boundary condition
switch regime

    case "traction-free"
        surfaceTraction = 0;
        % Plastic strain remains constant
        f(2*nbin+1) = ep_dot; 

    case "elastic-shell"
        sig_trial = params.A*(params.alph*u(end) - ep);
        surfaceTraction = -2*params.delta*sig_trial;
        % Plastic strain remains constant
        f(2*nbin+1) = ep_dot;

    case "plastic-shell"
        sig_trial = params.A*(params.alph*u(end) - ep);
        surfaceTraction = -2*params.delta*sign(sig_trial)*params.sig_y;

        % Plastic consistency condition
        u_surf_dot = wdot(nbin);
        f(2*nbin+1) = ep_dot - params.alph*u_surf_dot; 
        % Extra alpha term due to scaling

    otherwise
        error("Unknown boundary regime: %s", regime);
end

% Ghost value chosen so the interpolated surface stress equals
% surfaceTraction:
sigma_rr(end) = 2*surfaceTraction - sigma_rr(end-1);

%% Mechanical Equilibrium residual 
sigma_rr_edge = (sigma_rr(2:end)+sigma_rr(1:end-1))./2;

sigma_tt_edge = (sigma_tt(2:end)+sigma_tt(1:end-1))./2;

sigma_tt_edge(end+1) = sigma_tt_edge(nbin-1);   %nbin-1?


f(1:nbin) = diff(sigma_rr)./params.dr + (2./params.r_edge(2:end)).*(sigma_rr_edge - sigma_tt_edge);

sigma_kk = (sigma_rr(1:end-1) + 2.*sigma_tt); 


%% Stress-dependent chemical potential

J = [0;
    params.r_edge(2:end-1).^2.*(diff(c)./params.dr);
    0]; 
laplacian_c = diff(J)./(params.dr.*params.r_bin.^2); 


mu = mu_chem(c) - params.kappa*laplacian_c  - params.gamm.*sigma_kk./3;

%% Cahn-Hilliard transport residual
c_edge = (c(2:end) + c(1:end-1))./2; %interpolates for c on edges

Fr = zeros(nbin+1,1);
Fr(2:end-1) = -params.D0.*c_edge.*(1-c_edge).*(diff(mu)./params.dr);

% Boundary conditions
Fr(1) = 0;
Fr(end) = I;

cellVolume = (4*pi/3).*(params.r_edge(2:end).^3 - params.r_edge(1:end-1).^3);

div_flux = 4*pi*params.r_edge(2:end).^2.*Fr(2:end) ...
    - 4*pi*params.r_edge(1:end-1).^2.*Fr(1:end-1);

f(nbin+1:2*nbin) = cellVolume.*c_dot + div_flux;


end

