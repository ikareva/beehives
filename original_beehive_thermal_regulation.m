%% Bee Hive Thermoregulation, original model 
%
% Paper: Adaptive receptor expression and the emergence of disease as loss
% of signaling homeostasis by Irina Kareva
%
% Purpose
%   This script reproduces a temperature-only agent-based
%   thermoregulation model of a honeybee colony. Individual bees follow simple
%   threshold rules (too cold / too hot) and contribute heating or cooling.
%   The hive temperature evolves via a lumped thermal balance driven by (i) exchange
%   with the external environment and (ii) the mean net effect of bee behaviors.
%
% Provenance and attribution
%   The original model concept and NetLogo implementation were developed by Scott E. Page.
%   The implementation here is a MATLAB reproduction adapted from NetLogo code and from
%   code provided by Scott E. Page (personal communication, 2005). This MATLAB version
%   is intended as a transparent, minimal reproduction of the temperature dynamics.
%
% How to run
%   - Run this script in MATLAB. It produces a single figure with hive and external
%     temperatures versus simulation step.
%
% Notes / limitations
%   - The model treats the hive as a single thermal compartment (no spatial gradients).
%   - Agent behavior is represented only through threshold-triggered heating/cooling
%     contributions; movement and spatial clustering are not explicitly modeled.
%   - Threshold distributions and behavioral gains are illustrative parameters and can
%     be modified in the CONFIG section below.
%
% Reproducibility
%   - Random draws of individual thresholds are made reproducible via rng(seed).
%
% -------------------------------------------------------------------------

%% Reproducibility
rng(7);                                         % Fix the random seed for reproducible threshold draws

%% ---------- CONFIG ----------
N_BEES = 120;                                   % Total number of bees (agents)
FRACTION_OLD = 0.45;                            % Fraction of older, endothermic bees
DT = 0.20;                                      % Time step for Euler update (arbitrary units)
TOTAL_STEPS = 700;                              % Total number of update steps to simulate

HEAT_CAPACITY = 45.0;                           % Thermal inertia; larger => slower temperature response
K_ENV = 0.60;                                   % Coupling strength to external environment (exchange term)

ALPHA_HEAT = 26.0;                              % Heating gain for older bees when cold (scaled contribution)
ALPHA_FAN  = 22.0;                              % Cooling gain for older bees when hot (scaled contribution)

T_HIVE_INIT = 33.0;                             % Initial hive temperature (°C)

% External temperature schedule: piecewise-constant values over time
% Each row is [start_step, temperature_C], where start_step is interpreted as 0-based.
EXT_SCHEDULE = [ ...
    0   20.0;   ...                             % From step 0: 20°C
    80  28.0;   ...                             % From step 80: 28°C
    160 34.0;   ...                             % From step 160: 34°C
    240 35.0;   ...                             % From step 240: 35°C (warm peak)
    320 32.0;   ...                             % From step 320: 32°C
    400 26.0;   ...                             % From step 400: 26°C
    467 10.0;   ...                             % From step 467: 10°C (cold drop)
    560 12.0;   ...                             % From step 560: 12°C
    640 15.0];                                  % From step 640: 15°C

%% ---------- Bees & thresholds ----------
% Partition bees into younger (less endothermic) and older (more endothermic) roles.
n_old   = floor(N_BEES * FRACTION_OLD);         % Number of older bees (integer)
n_young = N_BEES - n_old;                       % Number of younger bees

% Each bee has a "too-warm" threshold and a "too-cold" threshold.
% Warm thresholds differ by age class; cold thresholds are derived by subtracting a gap.
%
% Warm thresholds (°C):
%   Young: Normal(mean=35, sd=0.7), clipped to [33.5, 37.0]
%   Old:   Normal(mean=36, sd=0.7), clipped to [34.0, 38.5]
warm_y  = min(max(35.0 + 0.7*randn(n_young,1), 33.5), 37.0);
warm_o  = min(max(36.0 + 0.7*randn(n_old,1),   34.0), 38.5);
warm_th = [warm_y; warm_o];                     % Concatenate into an N_BEES-by-1 vector

% Cold thresholds (°C):
%   cold_th = warm_th - gap, where gap ~ Uniform(1.2, 2.0)
% This enforces a neutral comfort band between cold_th and warm_th for each bee.
cold_th = warm_th - (1.2 + (2.0-1.2)*rand(N_BEES,1));

% Role encoding:
%   0 = young bee
%   1 = old bee
roles = [zeros(n_young,1); ones(n_old,1)];

%% ---------- External temperature helper ----------
% Returns the external temperature at a given integer step based on EXT_SCHEDULE.
% The schedule is interpreted as piecewise-constant with change points at start_step.
function Te = ext_temp(step, schedule)
    Te = schedule(1,2);                         % Default to the first temperature
    for k = 1:size(schedule,1)                  % Iterate through schedule rows in ascending order
        if step >= schedule(k,1)                % If step is beyond this change point,
            Te = schedule(k,2);                 % adopt this temperature.
        else
            break;                              % Stop once the appropriate interval is found
        end
    end
end

%% ---------- Simulation ----------
% State variable:
%   T_hive: lumped hive temperature (°C)
%
% At each step:
%   1) Determine which bees perceive "too cold" or "too hot" based on current T_hive.
%   2) Compute the magnitude of discomfort (gap) for those bees.
%   3) Convert discomfort to heating/cooling contributions depending on role.
%   4) Update T_hive via a lumped thermal balance (Euler step).

T_hive = T_HIVE_INIT;                           % Initialize hive temperature
temps  = zeros(TOTAL_STEPS,1);                  % Record hive temperature over time
exts   = zeros(TOTAL_STEPS,1);                  % Record external temperature over time

for t = 1:TOTAL_STEPS                           % Main simulation loop
    T_ext = ext_temp(t-1, EXT_SCHEDULE);        % External temperature (schedule uses 0-based step indexing)
    exts(t) = T_ext;                            % Log external temperature

    % Identify bees whose thresholds are violated at current hive temperature.
    too_cold = (T_hive < cold_th);              % Bees that perceive it as too cold
    too_hot  = (T_hive > warm_th);              % Bees that perceive it as too hot

    % Discomfort magnitudes (nonnegative); zero for bees within comfort band.
    cold_gap = (cold_th - T_hive) .* too_cold;  % If too cold: positive gap; else 0
    hot_gap  = (T_hive - warm_th) .* too_hot;   % If too hot: positive gap; else 0

    % Convert discomfort to heating/cooling contributions per bee.
    heat = zeros(N_BEES,1);                     % Heating contribution per bee (arbitrary units)
    cool = zeros(N_BEES,1);                     % Cooling contribution per bee (arbitrary units)

    % Older bees contribute more strongly to both heating and cooling.
    heat(roles==1) = ALPHA_HEAT * cold_gap(roles==1);
    cool(roles==1) = ALPHA_FAN  * hot_gap(roles==1);

    % Younger bees contribute weakly (scaled relative to older bees).
    heat(roles==0) = 0.4  * cold_gap(roles==0);
    cool(roles==0) = 0.15 * hot_gap(roles==0);

    % Net contribution per bee; positive warms, negative cools.
    net_bee_flux = heat - cool;

    % Lumped thermal balance:
    %   Environmental exchange: K_ENV*(T_ext - T_hive)
    %   Behavioral forcing:     mean(net_bee_flux)
    %   Thermal inertia:        HEAT_CAPACITY
    %
    % Euler update:
    %   T_hive(t+1) = T_hive(t) + [(K_ENV*(T_ext - T_hive) + mean(net_bee_flux)) * DT / HEAT_CAPACITY]
    dT = (K_ENV*(T_ext - T_hive) + mean(net_bee_flux)) * (DT / HEAT_CAPACITY);
    T_hive = T_hive + dT;                       % Advance hive temperature
    temps(t) = T_hive;                          % Record hive temperature
end

%% ---------- Plot ----------
figure;                                         % New figure
plot(temps, 'DisplayName','Hive'); hold on;     % Hive temperature
plot(exts,  'DisplayName','External');          % External temperature
xlabel('Step'); ylabel('Temperature (°C)');     % Axes labels
title('Hive vs External Temperature');          % Title
legend('show');                                 % Legend
grid on;                                        % Optional: improves readability for reproduction plots