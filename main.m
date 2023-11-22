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

% Total slot duration
total_tau = 10e-3;

% Number of blocks to reduce complexity (e.g. if N_RIS = 100 and N_blocks =
% 100/4 = 25, 4 elements together are optimized)
N_blocks = N/2;     
% be very careful in choosing the right divider respect to N_RIS

% Number of bits used to quantize the phase of the RIS
qtz_RIS = 1; 

% Load receive combining (CHECK THE CORRESPONDING N_AP)
load(['data/combining/w_candidates_' num2str(M) '_antennas.mat'], 'weights')
possible_angles = size(weights, 2);

%% Prepare to save simulation results
total_energy = zeros(num_setups, 1);
avg_delay = zeros(num_setups, 1);

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

    % Optmize control
    [current_total_energy, current_avg_delay] = ...
    rismec_control(N_slot, K, M, N, total_tau, N_blocks, qtz_RIS, weights, possible_angles, ...
        current_channel_ap_ris, current_channel_ap_ue, current_channel_ue_ris);

    % Save current data
    total_energy(rr) = current_total_energy;
    avg_delay(rr) = current_avg_delay;

    simulation_time = toc;
    disp(['elapsed ', num2str(simulation_time), ' seconds.']);

end

% Save results
file_name = ['data/results/N' num2str(N), '_M', num2str(M), '_K', num2str(K), '.mat'];
save(file_name)

