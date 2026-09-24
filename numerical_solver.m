function raw = numerical_solver(params)

nbin = params.nbin;

durations = params.durations;
currents = params.currents;


u0 = 0*ones(1,nbin); %excludes origin


c0 = 0.05.*ones(1,nbin); 
w0 = [u0, c0,0];
wdot0 = zeros(1,length(w0));



t_total = [];
sol_total = [];
t_plast = [];
t_elast = [];



%RHS_current = @(t,w,wdot) RHS_elastic(t,w,wdot,params);
%options = options_elast;
isPlastic = false;

t0 = 0;
I_prev = 0;


for i=1:length(durations)
   
    tf = t0 + durations(i);
    [t_seg, sol_seg, tseg_plast, tseg_elast, isPlastic] = solve_segment([t0, tf],w0,wdot0,params,currents(i),I_prev, isPlastic);
    t_total = [t_total;t_seg];
    sol_total = [sol_total;sol_seg];
    t_plast = [t_plast,tseg_plast]
    t_elast = [t_elast,tseg_elast]
    t0 = t_total(end)+1;
    w0 = sol_seg(end,:);
    %calculate wdot0
    wdot0 = (sol_total(end,:) - sol_total(end-1,:))/(t_total(end)-t_total(end-1));

    I_prev = currents(i);
    
end


sol_u = sol_total(:,1:nbin);
sol_c = sol_total(:,nbin+1:2*nbin);
sol_ept = sol_total(:,2*nbin+1);


I = zeros(length(t_total),1);
counter = 1;
for i = 1:length(t_total) 
    if t_total(i) > cumsum(durations(counter))&& counter < length(durations)
        counter = counter +1;
    end
    I(i) = currents(counter);
end

raw.sol_total = sol_total;
raw.t_total = t_total;
raw.I = I;
raw.t_plast = t_plast;
raw.t_elast = t_elast;
end