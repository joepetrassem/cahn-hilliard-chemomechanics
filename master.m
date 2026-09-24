
clc;clearvars;
tic
params = parameters();  
%params.shellModel = "traction-free";
% params.shellModel = "elastic";
 params.shellModel = "elastoplastic";


raw = numerical_solver(params); 

sol = postProcessing(raw, params);
figs = visualisation(sol,params);

sol_save = struct();
sol_save.t = sol.t_total;
sol_save.soc = sol.SOC;
sol_save.mu = sol.mu(:,end);
sol_save.c = sol.c(:,1:params.nbin/100:end);
sol_save.u = sol.u(:,1:params.nbin/100:end);
sol_save.params = params;
toc
