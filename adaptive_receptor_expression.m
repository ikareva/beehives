%% Receptor Homeostasis — Expression Dynamics with Live Visualization (No Keyboard Controls)
%
% Paper: Adaptive receptor expression and the emergence of disease as loss of signaling homeostasis
%
% Purpose
%   This script simulates a minimal, agent-based representation of receptor expression
%   homeostasis under a time-varying external ligand signal. A single lumped “cell signaling”
%   variable is driven by the external ligand and by aggregate receptor actions. Individual
%   receptors adjust a normalized surface expression proxy e_i in response to whether the
%   global signaling state is below or above each receptor’s comfort band.
%
% Conceptual mapping (high-level)
%   - External ligand/input:           L_ext(t)
%   - Lumped cell signaling state:     S_cell(t)
%   - Receptors (agents):              i = 1..N_RECEPTORS
%   - Surface expression proxy:        e_i(t) in [0,1]
%   - Regulation:
%       S_cell too low  -> upregulate (increase e_i)
%       S_cell too high -> internalize (decrease e_i)
%       optional shedding -> e_i set to 0 permanently
%       within comfort band -> relax toward baseline e0_i
%
% Outputs
%   - Panel A: receptor locations (static for display) colored by regulation state
%   - Panel B: external ligand trajectory (live)
%   - Panel C: counts in each regulation state + total surface expression Σe (live)
%   - Panel D: lumped signaling trajectory (live), optionally with a “normative” band
%
% Repository notes
%   - This script is self-contained (local functions at the end).
%   - Random elements are controlled via rng(seed) for reproducibility.
%   - The ligand schedule is generated programmatically via scenario presets.
%   - Keyboard controls (pause/stop/speed) are intentionally not included to
%     avoid conflicts with MATLAB interaction modes.
%
% -------------------------------------------------------------------------

%% Reproducibility
rng(7);  % Reproducible random draws (thresholds, baseline expression, noise)

%% =========================== USER SETTINGS ===========================
% Model size / display geometry
N_RECEPTORS   = 120;        % total receptors (agents)
FRACTION_REG  = 0.45;       % fraction of "regulatory" receptors (stronger effect)
CELL_RADIUS   = 1.0;        % radius for display-only cell boundary

% Threshold heterogeneity
HETEROGENEOUS_THRESHOLDS = true;    % reserved for future use (uniform thresholds option)

% Simulation horizon and rendering cadence
SIM_DT = 0.20;              % time step for Euler update (arbitrary units)
TOTAL_STEPS = 1200;         % number of simulation steps
PHYS_STEPS_PER_FRAME = 1;   % simulate this many steps per rendered frame (higher = faster)
DRAW_EVERY = 1;             % draw every Nth step (higher = fewer renders)
PAUSE_FACTOR = 0.00;        % optional slowdown factor per rendered frame (0 = no extra pause)

% Lumped signaling dynamics
SIGNAL_CAPACITY = 45.0;     % inertia; larger => slower response
K_LIGAND        = 0.60;     % coupling of signaling to ligand

% Aggregate action gains (analogue of “heating/cooling” strengths)
BETA_UP   = 26.0;           % gain applied when signaling is too low (upregulation pressure)
BETA_DOWN = 22.0;           % gain applied when signaling is too high (internalization pressure)

% Expression kinetics for each receptor (no spatial migration)
EXPR_RATE_UP   = 0.40;      % rate for increasing e when upregulating (scaled by SIM_DT)
EXPR_RATE_DOWN = 0.40;      % rate for decreasing e when internalizing (scaled by SIM_DT)
RELAX_RATE     = 0.08;      % relaxation rate toward baseline e0 when within comfort band
NOISE_LEVEL    = 0.04;      % additive noise magnitude for expression updates (0 = deterministic)

% Shedding options
SHEDDING_MODE = "conditional";   % "off" | "conditional" | "probabilistic"
SHED_MIN_E    = 0.02;            % shedding eligibility threshold on e (near-zero expression)
SHED_BIAS     = 0.8;             % additional signal above hi_thr that favors shedding
SHED_PROB     = 0.10;            % probability per step if eligible (probabilistic mode)

%% =================== EXTERNAL LIGAND: SCENARIO BUILDER ===================
% Scenario name (case-insensitive):
%   original_schedule
%   extreme_oscillation
%   mild_oscillation
%   consistently_low
%   consistently_high
%   slow_ramp
%   abrupt_jumps
%   chronic_increase
%   chronic_increase_with_intervention
SCENARIO = "chronic_increase";

% Scenario parameters (shared)
BASELINE      = 30;      % a.u., center value for oscillations
AMP_EXTREME   = 12;      % a.u., amplitude for extreme oscillation
AMP_MILD      = 4;       % a.u., amplitude for mild oscillation
PERIOD_STEPS  = 120;     % steps per full oscillation cycle
JITTER_SIGMA  = 0.0;     % a.u., add N(0, sigma^2) noise to ligand after scenario generation

% Consistently low/high presets
LOW_BASELINE  = 4;       % a.u.
HIGH_BASELINE = 42;      % a.u.
BASE_WOBBLE   = 2;       % a.u.

% Slow ramp preset
RAMP_UP_STEPS   = 200;   % steps to ramp up
RAMP_DOWN_STEPS = 200;   % steps to ramp down
LOW_LEVEL       = 16;    % start/end level of ramp
HIGH_LEVEL      = 38;    % top level of ramp
PLATEAU_STEPS   = 100;   % hold at high level

% Abrupt jumps preset
JUMP_LEVELS   = [18 34 22 36 20 30];  % sequence of plateaus
DWELL_STEPS   = 80;                   % steps per plateau

% Chronic increase presets
CI_LOWER = 16;          % a.u., starting/low bound
CI_UPPER = 44;          % a.u., max/high bound (pre-intervention ceiling)
CI_RATE  = 0.06;        % a.u. per step (pre-intervention slope)

% Chronic increase with intervention presets
INTV_STEP       = 320;           % step at which intervention occurs
INTV_MODE       = "absolute";    % "absolute" | "relative"
INTV_DROP_AU    = 40;            % absolute drop (a.u.) if INTV_MODE="absolute"
INTV_DROP_FRAC  = 0.40;          % fractional drop if INTV_MODE="relative" (0.4 = 40%)
INTV_FLOOR      = CI_LOWER;      % clamp after drop (do not go below this)
INTV_RESUME     = true;          % resume increasing after drop
CI_RATE_POST    = CI_RATE;       % post-intervention slope
CI_UPPER_POST   = CI_UPPER;      % post-intervention ceiling

% Build per-step ligand vector L_ext_vec for steps 0..TOTAL_STEPS-1
tSteps = (0:TOTAL_STEPS-1)';     % step index (column vector)
L_ext_vec = build_ligand_schedule( ...
    SCENARIO, tSteps, ...
    BASELINE, AMP_EXTREME, AMP_MILD, PERIOD_STEPS, ...
    LOW_BASELINE, HIGH_BASELINE, BASE_WOBBLE, ...
    RAMP_UP_STEPS, RAMP_DOWN_STEPS, LOW_LEVEL, HIGH_LEVEL, PLATEAU_STEPS, ...
    JUMP_LEVELS, DWELL_STEPS, ...
    CI_LOWER, CI_UPPER, CI_RATE, ...
    INTV_STEP, INTV_MODE, INTV_DROP_AU, INTV_DROP_FRAC, INTV_FLOOR, ...
    INTV_RESUME, CI_RATE_POST, CI_UPPER_POST);

% Optional global noise
if JITTER_SIGMA > 0
    L_ext_vec = L_ext_vec + JITTER_SIGMA * randn(TOTAL_STEPS,1);
end

% Optional rounding for display readability
L_ext_vec = round(L_ext_vec, 1);

% Schedule representation: [step, value] for piecewise-constant lookup
LIGAND_SCHEDULE = [tSteps, L_ext_vec];

%% =================== END EXTERNAL LIGAND: SCENARIO BUILDER ===================

%% Plot / figure formatting settings
TITLE_FONTSIZE = 14;
SUB_FONTSIZE   = 10;

% Display options
LEGEND_OUTSIDE   = true;          % horizontal legend below subsets plot
SHOW_NORM_BAND   = true;          % show a shaded “normative” band on signaling plot
NORM_BAND        = [32 36];       % a.u. (min, max)
NORM_BAND_COLOR  = [0.2 0.6 1.0]; % RGB
NORM_BAND_ALPHA  = 0.15;          % transparency (0..1)

%% =========================== INITIALIZATION ===========================
clip = @(x,a,b) min(max(x,a),b);  % clamp helper

% Agent classes: baseline vs regulatory
n_reg  = floor(N_RECEPTORS * FRACTION_REG);
n_base = N_RECEPTORS - n_reg;
classes = [zeros(n_base,1); ones(n_reg,1)];     % 0 = baseline, 1 = regulatory

% Comfort-band thresholds (per receptor)
hi_thr_base = clip(35.0 + 0.7*randn(n_base,1), 33.5, 37.0);
hi_thr_reg  = clip(36.0 + 0.7*randn(n_reg,1),  34.0, 38.5);
hi_thr = [hi_thr_base; hi_thr_reg];

% "Too-low" thresholds spaced below hi_thr by a random gap
lo_thr = hi_thr - (1.2 + (2.0-1.2)*rand(N_RECEPTORS,1));

% Baseline expression and initialization
e0_base = 0.35 + 0.10*rand(n_base,1);
e0_reg  = 0.55 + 0.10*rand(n_reg,1);
e0      = [e0_base; e0_reg];
e       = e0;

% Static receptor positions for display only (no migration)
r_disp_base = 0.35 + 0.15*rand(n_base,1);
r_disp_reg  = 0.60 + 0.20*rand(n_reg,1);
r_disp      = [r_disp_base; r_disp_reg];
theta       = 2*pi*rand(N_RECEPTORS,1);
[x_disp, y_disp] = pol2cart(theta, CELL_RADIUS * r_disp);

% Marker size encodes hi_thr
size_min = 40; size_max = 160;
sizes = size_min + (hi_thr - min(hi_thr)) ./ (max(hi_thr)-min(hi_thr)+eps) * (size_max - size_min);

% Lumped signaling state
S_cell = 33.0;

% Logs
S_log  = zeros(TOTAL_STEPS,1); % retained for potential downstream saving
L_log  = zeros(TOTAL_STEPS,1); 
counts_inband   = zeros(TOTAL_STEPS,1);
counts_up       = zeros(TOTAL_STEPS,1);
counts_internal = zeros(TOTAL_STEPS,1);
counts_shed     = zeros(TOTAL_STEPS,1);
total_expr      = zeros(TOTAL_STEPS,1);

% Permanent shedding mask
isShed = false(N_RECEPTORS,1);

%% =========================== FIGURE SETUP ===========================
tl = tiledlayout(2,2,"TileSpacing","compact","Padding","compact");

% Panel A: receptor display (static positions; colors change by state)
axCell = nexttile(tl,1);
axis(axCell,'equal'); hold(axCell,'on');
xlim(axCell,[-1.15 1.5]); ylim(axCell,[-1.25 1.05]); axis(axCell,'off');

text(axCell, 0.5, 1.12, 'Receptor Homeostasis (Expression Dynamics)', ...
    'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom', ...
    'FontSize',TITLE_FONTSIZE,'FontWeight','bold','Clipping','off');

text(axCell, 0.5, 1.06, 'In-band: blue  |  Up: green  |  Internalized: orange  |  Shed: gray', ...
    'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom', ...
    'FontSize',SUB_FONTSIZE,'Clipping','off');

% Cell boundary
tt = linspace(0,2*pi,400);
plot(axCell, CELL_RADIUS*cos(tt), CELL_RADIUS*sin(tt), 'k', 'LineWidth',2);

% Scatter object (colors updated during simulation)
hReceptors = scatter(axCell, x_disp, y_disp, sizes, 'filled', ...
    'MarkerEdgeColor','k','LineWidth',0.6);

% HUD beneath panel A
hudTxt = text(axCell, 0.5, -0.01, '', 'Units','normalized','Clipping','off', ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',12,'FontWeight','bold','Color',[0 0 0], ...
    'BackgroundColor',[0.95 0.95 0.95],'EdgeColor','k','Margin',6);

% Panel B: ligand (live)
axLig = nexttile(tl,2); hold(axLig,'on'); grid(axLig,'on');
xlim(axLig,[0 TOTAL_STEPS]);
yMinSched = min(LIGAND_SCHEDULE(:,2)); yMaxSched = max(LIGAND_SCHEDULE(:,2));
ylim(axLig,[yMinSched-5, yMaxSched+5]);
xlabel(axLig,'Step'); ylabel(axLig,'Ligand (a.u.)');
title(axLig,'External Ligand (live)');
hLigLine = animatedline(axLig,'LineWidth',1.5);

% Panel C: subsets + total expression
axSub = nexttile(tl,3); hold(axSub,'on'); grid(axSub,'on');
xlim(axSub,[0 TOTAL_STEPS]);
title(axSub,'Expression Subsets (live)');

yyaxis(axSub,'left');
ylabel(axSub,'Count (agents)');
hLnInBand = animatedline(axSub,'Color',[0 0.35 0.90], 'LineWidth',1.5, 'DisplayName','In-band');
hLnUp     = animatedline(axSub,'Color',[0 0.60 0.00], 'LineWidth',1.5, 'DisplayName','Upregulated');
hLnInt    = animatedline(axSub,'Color',[0.90 0.50 0.00], 'LineWidth',1.5, 'DisplayName','Internalized');
hLnShed   = animatedline(axSub,'Color',[0.50 0.50 0.50], 'LineWidth',1.5, 'DisplayName','Shed');

yyaxis(axSub,'right');
ylabel(axSub,'Total surface expression (Σe)');
hLnTotal  = animatedline(axSub,'LineStyle','--','LineWidth',1.5,'DisplayName','Total Σe');

yyaxis(axSub,'left');

if LEGEND_OUTSIDE
    legend(axSub,'Orientation','horizontal','Location','southoutside','Box','off');
else
    legend(axSub,'Location','northwest');
end

% Panel D: signaling (live)
axSig = nexttile(tl,4); hold(axSig,'on'); grid(axSig,'on');
xlim(axSig,[0 TOTAL_STEPS]);
ylim(axSig,[yMinSched-5, yMaxSched+5]);
xlabel(axSig,'Step'); ylabel(axSig,'Signaling (a.u.)');
title(axSig,'Cell Signaling (live)');
ylim(axSig,[20 50]);

if SHOW_NORM_BAND
    xBand = [0, TOTAL_STEPS, TOTAL_STEPS, 0];
    yBand = [NORM_BAND(1), NORM_BAND(1), NORM_BAND(2), NORM_BAND(2)];
    hNorm = patch('Parent',axSig, 'XData',xBand, 'YData',yBand, ...
        'FaceColor',NORM_BAND_COLOR, 'FaceAlpha',NORM_BAND_ALPHA, ...
        'EdgeColor','none', 'HandleVisibility','off');
    uistack(hNorm,'bottom');
end

hSigLine = animatedline(axSig,'LineWidth',1.5);

% Panel labels
LABEL_FONTSIZE = 16;
LABEL_FONTWEIGHT = 'bold';
LABEL_COLOR = 'k';

text(axCell, 0.02, 0.02, 'A', 'Units','normalized', ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
    'FontSize',LABEL_FONTSIZE, 'FontWeight',LABEL_FONTWEIGHT, ...
    'Color',LABEL_COLOR, 'Clipping','off');

text(axLig, 0.02, 0.02, 'B', 'Units','normalized', ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
    'FontSize',LABEL_FONTSIZE, 'FontWeight',LABEL_FONTWEIGHT, ...
    'Color',LABEL_COLOR);

text(axSub, 0.02, 0.02, 'C', 'Units','normalized', ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
    'FontSize',LABEL_FONTSIZE, 'FontWeight',LABEL_FONTWEIGHT, ...
    'Color',LABEL_COLOR);

text(axSig, 0.02, 0.02, 'D', 'Units','normalized', ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
    'FontSize',LABEL_FONTSIZE, 'FontWeight',LABEL_FONTWEIGHT, ...
    'Color',LABEL_COLOR);

%% =========================== SIMULATION LOOP ===========================
state_draw = zeros(N_RECEPTORS,1);
L_ext_draw = LIGAND_SCHEDULE(1,2);

for t = 1:TOTAL_STEPS

    % Multiple physics steps per rendered frame
    for k = 1:PHYS_STEPS_PER_FRAME
        step_idx = min(t + k - 1, TOTAL_STEPS);

        % External ligand at this step
        L_ext = ligand_at(step_idx-1, LIGAND_SCHEDULE);   % schedule uses 0-based step indexing
        L_ext_draw = L_ext;

        % Comfort-band checks (based on the current global signaling state)
        too_low  = (S_cell < lo_thr);
        too_high = (S_cell > hi_thr);

        low_gap  = (lo_thr - S_cell) .* double(too_low);
        high_gap = (S_cell - hi_thr) .* double(too_high);

        % Aggregate regulatory “action” (used only as a forcing term on S_cell)
        up   = zeros(N_RECEPTORS,1);
        down = zeros(N_RECEPTORS,1);

        up(classes==1)   = BETA_UP   * low_gap(classes==1);
        down(classes==1) = BETA_DOWN * high_gap(classes==1);

        up(classes==0)   = 0.4  * low_gap(classes==0);
        down(classes==0) = 0.15 * high_gap(classes==0);

        net_action = up - down;

        % Lumped signaling dynamics (Euler step)
        dS = (K_LIGAND*(L_ext - S_cell) + mean(net_action)) * (SIM_DT / SIGNAL_CAPACITY);
        S_cell = S_cell + dS;

        % Expression dynamics (per receptor)
        state = zeros(N_RECEPTORS,1);        % 0 in-band, +1 upregulated, -1 internalized, -2 shed
        alive = ~isShed;

        % Upregulate: increase e when signaling is too low
        idx = alive & too_low;
        if any(idx)
            e(idx) = clip(e(idx) + EXPR_RATE_UP*SIM_DT + NOISE_LEVEL*randn(nnz(idx),1), 0, 1);
            state(idx) = 1;
        end

        % Internalize: decrease e when signaling is too high
        idx = alive & too_high;
        if any(idx)
            e(idx) = e(idx) - EXPR_RATE_DOWN*SIM_DT + NOISE_LEVEL*randn(nnz(idx),1);
            e(idx) = clip(e(idx), 0, 1);
            state(idx) = -1;
        end

        % In-band: relax toward baseline expression e0
        idx = alive & ~(too_low | too_high);
        if any(idx)
            drift = RELAX_RATE * SIM_DT * sign(e0(idx) - e(idx));
            e(idx) = clip(e(idx) + drift + NOISE_LEVEL*randn(nnz(idx),1), 0, 1);
        end

        % Shedding logic (optional)
        switch lower(SHEDDING_MODE)
            case 'off'
                shed_now = false(N_RECEPTORS,1);

            case 'conditional'
                eligible = alive & (e <= SHED_MIN_E) & (S_cell > hi_thr + SHED_BIAS);
                shed_now = eligible;

            case 'probabilistic'
                eligible = alive & ((e <= SHED_MIN_E) | (S_cell > hi_thr + SHED_BIAS));
                shed_now = eligible & (rand(N_RECEPTORS,1) < SHED_PROB);

            otherwise
                error('Unknown SHEDDING_MODE: %s', SHEDDING_MODE);
        end

        if any(shed_now)
            isShed = isShed | shed_now;
            e(shed_now) = 0;
        end

        state(isShed) = -2;

        % Logs for subset plot and HUD
        counts_inband(step_idx)   = sum(state==0);
        counts_up(step_idx)       = sum(state==1);
        counts_internal(step_idx) = sum(state==-1);
        counts_shed(step_idx)     = sum(isShed);
        total_expr(step_idx)      = sum(e);

        % Retained logs (not strictly required for plotting)
        S_log(step_idx) = S_cell; 
        L_log(step_idx) = L_ext;  
    end

    % Render (throttled)
    if mod(t, DRAW_EVERY) == 0

        % Colors by state: in-band blue; up green; internalized orange; shed gray
        C = repmat([0.00 0.35 0.90], N_RECEPTORS, 1);
        idx = (state==1);   if any(idx), C(idx,:) = repmat([0.00 0.65 0.00], nnz(idx), 1); end
        idx = (state==-1);  if any(idx), C(idx,:) = repmat([0.95 0.55 0.00], nnz(idx), 1); end
        idx = (state==-2);  if any(idx), C(idx,:) = repmat([0.60 0.60 0.60], nnz(idx), 1); end
        set(hReceptors, 'CData', C);

        % Heads-up display
        hudStr = sprintf('Signaling: %.1f a.u.    Ligand: %.1f a.u.    Total surface expr: %.1f    Shedding: %s', ...
            S_cell, L_ext_draw, sum(e), upper(SHEDDING_MODE));
        set(hudTxt,'String', hudStr);

        % Live plots
        addpoints(hLigLine, t, L_ext_draw);
        addpoints(hSigLine, t, S_cell);

        yyaxis(axSub,'left');
        addpoints(hLnInBand, t, counts_inband(t));
        addpoints(hLnUp,     t, counts_up(t));
        addpoints(hLnInt,    t, counts_internal(t));
        addpoints(hLnShed,   t, counts_shed(t));

        yyaxis(axSub,'right');
        addpoints(hLnTotal,  t, total_expr(t));

        drawnow limitrate nocallbacks;

        if PAUSE_FACTOR > 0
            pause(PAUSE_FACTOR * SIM_DT);
        end
    end
end

%% ========================= LOCAL FUNCTIONS ===========================
function L = ligand_at(step, schedule)
% ligand_at  Piecewise-constant lookup from an [step, value] schedule.
%           step is interpreted as 0-based.
    L = schedule(1,2);
    for k = 1:size(schedule,1)
        if step >= schedule(k,1)
            L = schedule(k,2);
        else
            break;
        end
    end
end

function L_ext_vec = build_ligand_schedule( ...
    SCENARIO, t, ...
    BASELINE, AMP_EXTREME, AMP_MILD, PERIOD_STEPS, ...
    LOW_BASELINE, HIGH_BASELINE, BASE_WOBBLE, ...
    RAMP_UP_STEPS, RAMP_DOWN_STEPS, LOW_LEVEL, HIGH_LEVEL, PLATEAU_STEPS, ...
    JUMP_LEVELS, DWELL_STEPS, ...
    CI_LOWER, CI_UPPER, CI_RATE, ...
    INTV_STEP, INTV_MODE, INTV_DROP_AU, INTV_DROP_FRAC, INTV_FLOOR, ...
    INTV_RESUME, CI_RATE_POST, CI_UPPER_POST)

    TOTAL_STEPS = numel(t);
    L_ext_vec = zeros(TOTAL_STEPS,1);

    switch lower(string(SCENARIO))

        case "original_schedule"
            ORIG_SCHEDULE = [ ...
                0   20.0;
                80  28.0;
                160 34.0;
                240 35.0;
                320 32.0;
                400 26.0;
                467 10.0;
                560 12.0;
                640 15.0];

            L_ext_vec = ORIG_SCHEDULE(1,2) * ones(TOTAL_STEPS,1);
            for k = 1:size(ORIG_SCHEDULE,1)
                step_start = ORIG_SCHEDULE(k,1);
                val = ORIG_SCHEDULE(k,2);
                L_ext_vec(t >= step_start) = val;
            end

        case "extreme_oscillation"
            s = sign(sin(2*pi*t/PERIOD_STEPS));
            L_ext_vec = BASELINE + AMP_EXTREME * s;

        case "mild_oscillation"
            L_ext_vec = BASELINE + AMP_MILD * sin(2*pi*t/PERIOD_STEPS);

        case "consistently_low"
            L_ext_vec = LOW_BASELINE + BASE_WOBBLE * sin(2*pi*t/(2*PERIOD_STEPS));

        case "consistently_high"
            L_ext_vec = HIGH_BASELINE + BASE_WOBBLE * sin(2*pi*t/(2*PERIOD_STEPS));

        case "slow_ramp"
            L_ext_vec = LOW_LEVEL * ones(TOTAL_STEPS,1);
            idx = 1;

            n = min(RAMP_UP_STEPS, TOTAL_STEPS-idx+1);
            if n>0
                L_ext_vec(idx:idx+n-1) = linspace(LOW_LEVEL, HIGH_LEVEL, n);
                idx = idx + n;
            end

            n = min(PLATEAU_STEPS, TOTAL_STEPS-idx+1);
            if n>0
                L_ext_vec(idx:idx+n-1) = HIGH_LEVEL;
                idx = idx + n;
            end

            n = min(RAMP_DOWN_STEPS, TOTAL_STEPS-idx+1);
            if n>0
                L_ext_vec(idx:idx+n-1) = linspace(HIGH_LEVEL, LOW_LEVEL, n);
                idx = idx + n;
            end

            if idx <= TOTAL_STEPS
                L_ext_vec(idx:end) = LOW_LEVEL;
            end

        case "abrupt_jumps"
            idx = 1; lev = 1;
            while idx <= TOTAL_STEPS
                n = min(DWELL_STEPS, TOTAL_STEPS-idx+1);
                L_ext_vec(idx:idx+n-1) = JUMP_LEVELS(lev);
                idx = idx + n;
                lev = lev + 1;
                if lev > numel(JUMP_LEVELS)
                    lev = 1;
                end
            end

        case "chronic_increase"
            L_ext_vec = CI_LOWER + CI_RATE * t;
            L_ext_vec = min(L_ext_vec, CI_UPPER);

        case "chronic_increase_with_intervention"
            L_pre = CI_LOWER + CI_RATE * t;
            L_pre = min(L_pre, CI_UPPER);

            t0 = max(0, min(TOTAL_STEPS-1, INTV_STEP));
            L0_pre = L_pre(t0+1);

            switch lower(string(INTV_MODE))
                case "absolute"
                    L0_drop = max(INTV_FLOOR, L0_pre - INTV_DROP_AU);
                case "relative"
                    L0_drop = max(INTV_FLOOR, L0_pre * (1 - INTV_DROP_FRAC));
                otherwise
                    error('Unknown INTV_MODE: %s', INTV_MODE);
            end

            if INTV_RESUME
                L_ext_vec = L_pre;
                post = (t >= t0);
                L_ext_vec(post) = L0_drop + CI_RATE_POST * (t(post) - t0);
                L_ext_vec(post) = min(L_ext_vec(post), CI_UPPER_POST);
            else
                L_ext_vec = L_pre;
                L_ext_vec(t >= t0) = L0_drop;
            end

        otherwise
            error('Unknown SCENARIO: %s', SCENARIO);
    end
end