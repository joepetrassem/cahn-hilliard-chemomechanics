function [t_total,sol_total,t_plast,t_elast, isPlastic] = solve_segment(tspan,w0, wdot0,params,I,I_prev, isPlastic)

nbin = params.nbin;

t0 = tspan(1);
tf = tspan(2);


t_total = [];
sol_total = [];
t_plast = [];
t_elast = [];


options = odeset("RelTol",1e-7,"AbsTol",1e-9);

if isPlastic && I*I_prev<=0
    isPlastic = ~isPlastic;   %makes the system elastic if the current direction changes and was plastic
end


switch params.shellModel
    case "traction-free"
        regime = "traction-free";
        options.Events = [];

    case "elastic"
        regime = "elastic-shell";
        options.Events = [];  % Never switch to plastic

    case "elastoplastic"
        if isPlastic
            regime = "plastic-shell";
            options.Events = [];
        else
            regime = "elastic-shell";
            options.Events = ...
            @(t,w,wdot) elastEvent(t,w,wdot,params);
        end
        
    otherwise
        error("Unknown shell model: %s", params.shellModel);

end






RHS_current = @(t,w,wdot) RHS_particle(t,w,wdot,params, I,regime);


while t0 < tf
    [w0, wdot0] = decic(RHS_current,t0,w0,[ones(1,nbin), zeros(1,nbin+1)],wdot0,[]);

    eventFcn = odeget(options,'Events');

    if isempty(eventFcn) %no elastoplastic shell case
        [t_seg, sol_seg] = ode15i(RHS_current, [t0 tf], w0, wdot0, options);
    
        te = [];    
    else                 %elastoplastic shell case
        [t_seg, sol_seg, te, we, ie] = ode15i(RHS_current, [t0 tf], w0, wdot0, options);
    end
    
    
    t_total = [t_total; t_seg];
    sol_total = [sol_total; sol_seg];
    
    t0 = t_seg(end) + 1;
    w0 = sol_seg(end,:);
    wdot0 = (sol_total(end,:) - sol_total(end-1,:))/(t_total(end)-t_total(end-1));


    if ~isempty(te)  %Triggered elastEvent, i.e. yielded, went plastic
        isPlastic = true;
        RHS_current = @(t,w,wdot) RHS_particle(t,w,wdot,params, I, 'plastic-shell');
        t_plast = [t_plast,te(1)];

    elseif isPlastic && t_total(end) == tf
        t_elast = [t_elast, t_total(end)];

    end

end



end