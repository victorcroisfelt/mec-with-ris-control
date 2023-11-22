clear all; clc; rng('default');

%% Load Channel and Simulation Parameters 
file_name = 'data/channel/N64_M8_K2.mat';
load(file_name)

disp('----- Main Parameters -----')
disp(['N_slot = ', num2str(N_slot)])
disp(['K = ', num2str(K)])
disp(['N = ', num2str(N)])
disp(['M = ', num2str(N)])
disp(['num_setups = ', num2str(num_setups)])
disp('---------------------------')

% Number of blocks to reduce complexity (e.g. if N_RIS = 100 and N_blocks =
% 100/4 = 25, 4 elements together are optimized)
N_blocks = N/2;     
% be very careful in choosing the right divider respect to N_RIS

% Number of bits used to quantize the phase of the RIS
qtz_RIS = 1; 

% Load receive combining (CHECK THE CORRESPONDING N_AP)
load(['data/combining/w_candidates_' num2str(M) '_antennas.mat'], 'weights')
possible_angles = size(weights, 2);


%% Computing Time

% Define total slot duration
tau = 400e-3;

% Define TTI time
constant = 4;
tti_time = constant * 1/14 * 10e-3;

% Define the RIS switching time
tau_ris = 0; % 1/14 * 10e-3

% Compute signaling time
tau_sig = 2 * tau_ris + 5 * tti_time; 

% Compute channel estimation overhead
tau_ce = (tau_ris + tti_time) * (N_blocks + 1);

% Compute mu constant
mu = 5 * K * (6 * 1 * (1 + 5 * N^2) + 5 * (1 + 4 * N^2)) + N;

% Compute RA time 
tau_ra = (4 + 1) * N_blocks * mu;
tau_ra = tau_ra / 4.5e9;

% Compute overhead time
tau_overhead = tau_sig + tau_ce + tau_ra;

% Percentage of slot duration dedicated to offloading
perc = 1 - (tau_overhead / tau);

%% Simulation Parameters

% Define error type
% 1: INI-U
% 2: INI-R
% 3: SET-U
% 4: SET-R
error_type = 'ini-u';
%error_type = 'ini-r';
%error_type = 'set-u';
%error_type = 'set-r';

% Probability vector to lose of the packets:
probas = [1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99];

% Get length of the probability vector
num_points = length(probas);

%% Prepare to save simulation results
total_energy = zeros(num_setups, num_points);
avg_delay = zeros(num_setups, num_points);

%% Simulation
disp('------- Simulation --------')

% Go through all setups
for rr = 1:num_setups

    disp(['setup = ', num2str(rr), ])
    tic;

    % Get current channels
    current_channel_ap_ris = squeeze(overall_channel_ap_ris(rr, :, :, :, :));
    current_channel_ap_ue = squeeze(overall_channel_ap_ue(rr, :, :, :));
    current_channel_ue_ris = squeeze(overall_channel_ue_ris(rr, :, :, :));

    % Go through all probability points
    for nn = 1:num_points
       
        disp(['    point = ', num2str(nn), ])

        % Current probability
        proba = probas(nn);

        % Optmize control
        [current_total_energy, current_avg_delay] = ...
        rismec_control(N_slot, K, M, N, tau, perc, N_blocks, qtz_RIS, weights, possible_angles, ...
            current_channel_ap_ris, current_channel_ap_ue, current_channel_ue_ris, proba, error_type);

        % Save current data
        total_energy(rr, nn) = current_total_energy;
        avg_delay(rr, nn) = current_avg_delay;

    end

    simulation_time = toc;
    disp(['elapsed ', num2str(simulation_time), ' seconds.']);

end

% Save results
file_name = ['data/results/', error_type, '/N', num2str(N), '_M', num2str(M), '_K', num2str(K), '.mat'];
save(file_name)

