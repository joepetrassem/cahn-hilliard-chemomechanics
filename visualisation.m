function figs = visualisation(sol, params)
%VISUALISE_SOLUTION Plot the main outputs of the particle model.
%
%   figs = visualise_solution(sol, params)
%
%   Inputs
%   ------
%   sol    Solution structure returned by numerical_solver and, where
%          applicable, postprocess_solution.
%   params Model parameter structure.
%
%   Output
%   ------
%   figs   Structure containing the generated figure handles.
%
%   Expected shell modes:
%       params.shellModel = "traction-free"
%       params.shellModel = "elastic"
%       params.shellModel = "elastoplastic"
%



mode = params.shellModel;
t = sol.t_total;
nt = numel(t);
nProfiles = min(6, nt);
profileIdx = unique(round(linspace(1, nt, nProfiles)));
profileLabels = compose('t = %.3g', t(profileIdx));

soc = sol.SOC;
I = sol.I;
ep = sol.ep;

figs = struct();

%% Overall simulation summary
figs.summary = figure( ...
    'Name', 'Simulation summary', ...
    'Color', 'w');

summaryLayout = tiledlayout(figs.summary, 2, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(summaryLayout);
plot(ax, t, soc, 'LineWidth', 1.5);
xlabel(ax, 'Time');
ylabel(ax, 'Average concentration');
format_axes(ax);

ax = nexttile(summaryLayout);
plot(ax, t, sol.c(:,end), 'LineWidth', 1.5);
xlabel(ax, 'Time');
ylabel(ax, 'Surface concentration');
format_axes(ax);

ax = nexttile(summaryLayout);
plot(ax, t, sol.u(:,end), 'LineWidth', 1.5);
xlabel(ax, 'Time');
ylabel(ax, 'Surface displacement');
format_axes(ax);

ax = nexttile(summaryLayout);
if ~isempty(I)
    plot(ax, t, I(:), 'LineWidth', 1.5);
    yline(ax, 0, ':k');
    xlabel(ax, 'Time');
    ylabel(ax, 'Applied flux/current');
    currentMax = max(abs(I));
    if currentMax > 0
        ylim(ax, 1.1.*[-currentMax, currentMax]);
    end
    format_axes(ax);
else
    unavailable_axis(ax, 'Current history not available');
end

ax = nexttile(summaryLayout);
if isfield(sol, 'mu') && ~isempty(sol.mu)
    plot(ax, t, sol.mu(:,end), 'LineWidth', 1.5);
    xlabel(ax, 'Time');
    ylabel(ax, 'Surface chemical potential');
    format_axes(ax);
elseif ~isempty(ep) && mode ~= "traction-free"
    plot(ax, t, ep(:), 'LineWidth', 1.5);
    xlabel(ax, 'Time');
    ylabel(ax, 'Shell plastic strain');
    format_axes(ax);
else
    unavailable_axis(ax, 'Chemical potential not available');
end


ax = nexttile(summaryLayout);
if isfield(sol, 'voltage') && ~isempty(sol.voltage)
    plot(ax, soc, sol.voltage(:), 'LineWidth', 1.5);
    xlabel(ax, 'Average concentration');
    ylabel(ax, 'Model voltage');
    format_axes(ax);
else
    unavailable_axis(ax, 'Voltage not available');
end


sgtitle(summaryLayout, sprintf('Simulation summary: %s', ...
    strrep(char(mode), '-', ' ')));

%% Radial profiles
figs.profiles = figure( ...
    'Name', 'Radial profiles', ...
    'Color', 'w');

profileLayout = tiledlayout(figs.profiles, 2, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(profileLayout);
plot_profile_set(ax, params.r_bin, sol.c, profileIdx);
xlabel(ax, 'Radius');
ylabel(ax, 'Concentration');
format_axes(ax);
legend(ax, profileLabels, 'Location', 'best');

ax = nexttile(profileLayout);
uWithOrigin = [zeros(nt,1), sol.u];
plot_profile_set(ax, params.r_edge, uWithOrigin, profileIdx);
xlabel(ax, 'Radius');
ylabel(ax, 'Radial displacement');
format_axes(ax);

ax = nexttile(profileLayout);
if isfield(sol, 'mu') && ~isempty(sol.mu)
    plot_profile_set(ax, params.r_bin, sol.mu, profileIdx);
    xlabel(ax, 'Radius');
    ylabel(ax, 'Chemical potential');
    format_axes(ax);
else
    unavailable_axis(ax, 'Chemical potential not available');
end

ax = nexttile(profileLayout);
if isfield(sol, 'sigma_rr') && ~isempty(sol.sigma_rr)
    plot_profile_set(ax, params.r_bin, sol.sigma_rr, profileIdx);
    xlabel(ax, 'Radius');
    ylabel(ax, 'Radial stress');
    format_axes(ax);
else
    unavailable_axis(ax, 'Radial stress not available');
end

ax = nexttile(profileLayout);
if isfield(sol, 'sigma_tt') && ~isempty(sol.sigma_tt)
    plot_profile_set(ax, params.r_bin, sol.sigma_tt, profileIdx);
    xlabel(ax, 'Radius');
    ylabel(ax, 'Hoop stress');
    format_axes(ax);
else
    unavailable_axis(ax, 'Hoop stress not available');
end

ax = nexttile(profileLayout);
if isfield(sol, 'sigma_kk') && ~isempty(sol.sigma_kk)
    plot_profile_set(ax, params.r_bin, sol.sigma_kk, profileIdx);
    xlabel(ax, 'Radius');
    ylabel(ax, 'Stress trace');
    format_axes(ax);
else
    unavailable_axis(ax, 'Stress trace not available');
end

sgtitle(profileLayout, sprintf('Radial fields: %s', ...
    strrep(char(mode), '-', ' ')));

%% Shell-specific quantities
if mode ~= "traction-free"
    if isempty(ep)
        warning('visualise_solution:MissingShellState', ...
            'Shell mode is active, but sol.ep or sol.ept is unavailable.');
        ep = zeros(nt,1);
    else
        ep = ep(:);
    end

    trialStress = params.A .* ...
        (params.alph .* sol.u(:,end) - ep);

    if isfield(sol, 'surfaceTraction') && ~isempty(sol.surfaceTraction)
        surfaceTraction = sol.surfaceTraction(:);
    elseif mode == "elastic"
        surfaceTraction = -2 .* params.delta .* trialStress;
    else
        shellStress = max(min(trialStress, params.sig_y), -params.sig_y);
        surfaceTraction = -2 .* params.delta .* shellStress;
    end

    figs.shell = figure( ...
        'Name', 'Shell response', ...
        'Color', 'w');

    shellLayout = tiledlayout(figs.shell, 3, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(shellLayout);
    plot(ax1, t, surfaceTraction, 'LineWidth', 1.5);
    ylabel(ax1, 'Surface traction');
    format_axes(ax1);

    ax2 = nexttile(shellLayout);
    plot(ax2, t, trialStress, 'LineWidth', 1.5, ...
        'DisplayName', 'Trial stress');
    hold(ax2, 'on');
    yline(ax2, params.sig_y, '--', 'Yield stress', ...
        'HandleVisibility', 'off');
    yline(ax2, -params.sig_y, '--', ...
        'HandleVisibility', 'off');
    hold(ax2, 'off');
    ylabel(ax2, 'Shell trial stress');
    format_axes(ax2);

    ax3 = nexttile(shellLayout);
    plot(ax3, t, ep, 'LineWidth', 1.5);
    xlabel(ax3, 'Time');
    ylabel(ax3, 'Shell plastic strain');
    format_axes(ax3);

    if mode == "elastoplastic"
        add_transition_lines(ax1, sol);
        add_transition_lines(ax2, sol);
        add_transition_lines(ax3, sol);
    end

    sgtitle(shellLayout, sprintf('Shell response: %s', ...
        strrep(char(mode), '-', ' ')));
end

drawnow;

end

function validate_inputs(sol, params)
requiredSol = {'u','c'};
requiredParams = {'r_bin','r_edge','R'};

for k = 1:numel(requiredSol)
    if ~isfield(sol, requiredSol{k})
        error('visualise_solution:MissingSolutionField', ...
            'Missing solution field: sol.%s', requiredSol{k});
    end
end

for k = 1:numel(requiredParams)
    if ~isfield(params, requiredParams{k})
        error('visualise_solution:MissingParameterField', ...
            'Missing parameter field: params.%s', requiredParams{k});
    end
end

if size(sol.u,1) ~= size(sol.c,1)
    error('visualise_solution:InconsistentTimeDimension', ...
        'sol.u and sol.c must contain the same number of time points.');
end
end

function t = get_time(sol)
if isfield(sol, 't')
    t = sol.t(:);
elseif isfield(sol, 't_total')
    t = sol.t_total(:);
else
    error('visualise_solution:MissingTime', ...
        'The solution must contain sol.t or sol.t_total.');
end

if numel(t) ~= size(sol.c,1)
    error('visualise_solution:InconsistentTime', ...
        'The time vector length must match the number of solution rows.');
end
end

function mode = get_shell_mode(params)
if isfield(params, 'shellModel')
    mode = lower(string(params.shellModel));
elseif isfield(params, 'shellEnabled')
    if params.shellEnabled
        mode = "elastoplastic";
    else
        mode = "traction-free";
    end
else
    mode = "traction-free";
end

switch mode
    case {"traction-free", "none", "no-shell"}
        mode = "traction-free";
    case {"elastic", "elastic-shell"}
        mode = "elastic";
    case {"elastoplastic", "elastoplastic-shell", "plastic"}
        mode = "elastoplastic";
    otherwise
        error('visualise_solution:UnknownShellMode', ...
            'Unknown shell mode: %s', mode);
end
end

function plot_profile_set(ax, radius, field, idx)
plot(ax, radius(:), field(idx,:).', 'LineWidth', 1.2);
end

function format_axes(ax)
grid(ax, 'on');
box(ax, 'on');
ax.FontSize = 11;
end

function unavailable_axis(ax, message)
axis(ax, 'off');
text(ax, 0.5, 0.5, message, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Units', 'normalized', ...
    'Color', [0.35 0.35 0.35]);
end

function add_transition_lines(ax, sol)
hold(ax, 'on');

if isfield(sol, 't_plast') && ~isempty(sol.t_plast)
    for k = 1:numel(sol.t_plast)
        xline(ax, sol.t_plast(k), '--r', 'Yield', ...
            'HandleVisibility', 'off');
    end
end

if isfield(sol, 't_elast') && ~isempty(sol.t_elast)
    for k = 1:numel(sol.t_elast)
        xline(ax, sol.t_elast(k), ':b', 'Elastic', ...
            'HandleVisibility', 'off');
    end
end

hold(ax, 'off');
end