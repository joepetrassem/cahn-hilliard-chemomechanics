function params = parameters()

%% Geometry and mesh

params.R = 1;
params.nbin = 800;

params.dr = params.R/params.nbin;

params.r_edge = linspace( ...
    0, params.R, params.nbin + 1)';

params.r_bin = linspace( ...
    params.dr/2, ...
    params.R - params.dr/2, ...
    params.nbin)';


%% Material properties

params.radius = 75e-9;
params.area = 4*pi*params.radius^2;

params.nu_si = 0.25; %Poisson silicon/particle
params.E_si = 60e9;  %Youngs modulus silicon/particle
params.E_bar_si = params.E_si/(1 + params.nu_si);

params.Far = 96500;
params.R_gT = 8.314*300;

params.D = 1e-3;

params.V_c = 8.2e-6;
params.c_m = 2.7e5;
params.c_ref = 0.05;

%% Dimensionless groups

params.nu = params.nu_si/(1 - 2*params.nu_si);

params.alph = 0.33;  %1/3 V_c*c_m
params.gamm = 2;     %V_c*E_bar_si/(R_gT);
params.D0 = 1e-3;    %D*area*c_m/(I0*radius); 
params.kappa = 1e-4; 


%% Shell model

params.shellModel = 'traction-free';
params.A = 1;            % E*(E+3*E_b)/(E+E_b); ratio of particle and shell stiffness        
params.delta = 0.1;
params.sig_y = 0.1;

params.E_sei = 2e9;
params.nu_sei = 0.3;
params.E_bar_sei = params.E_sei/(1 + params.nu_sei);
params.E_b_sei = params.E_bar_sei*params.nu_sei ...
    /(1 - 2*params.nu_sei);

params.E = params.E_bar_sei/params.E_bar_si;
params.E_b = params.E_b_sei/params.E_bar_si;



%% Initial and operating conditions
params.I0 = 1e-4;

params.durations = [150e8];%, 150e8];
params.currents =[-0.002e-8];%, 0.002e-8];

end