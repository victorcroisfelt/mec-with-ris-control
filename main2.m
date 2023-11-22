clear all; clc; rand('state',0); randn('state',0)

%% Load Channel and Simulation Parameters 
file_name = 'data/channel/K2.mat';
load(file_name)

disp('----- Main Parameters -----')
disp(['N_slot = ', num2str(N_slot)])
disp(['K = ', num2str(K)])
disp(['N_RIS = ', num2str(N_RIS)])
disp(['N_AP = ', num2str(Nr)])
disp(['num_setups = ', num2str(num_setups)])
disp('---------------------------')

% Total slot duration
total_delta = 10e-3;

% Number of blocks to reduce complexity (e.g. if N_RIS = 100 and N_blocks =
% 100/4 = 25, 4 elements together are optimized)
N_blocks = N_RIS/2;     
% be very careful in choosing the right divider respect to N_RIS

% Number of bits used to quantize the phase of the RIS
qtz_RIS = 1; 

% Load receive combining (CHECK THE CORRESPONDING N_AP)
load('data/combining/w_candidates_32_antennas.mat', 'weights')
possible_angles = size(weights, 2);

proba = [0.5, 0, 0, 0];


%% Prepare to save simulation results
total_energy = zeros(num_setups);
avg_delay = zeros(num_setups);

%% Simulation
disp('')
disp('------- Simulation --------')

% Go through all setups
for rr = 1:num_setups
    
    disp(['setup = ', num2str(rr)])
    tic;

    % Get current channels
    current_channel_ap_ris = squeeze(overall_channel_ap_ris(rr, :, :, :, :));
    current_channel_ap_ue = squeeze(overall_channel_ap_ue(rr, :, :, :));
    current_channel_ue_ris = squeeze(overall_channel_ue_ris(rr, :, :, :));

    % Optmize control
    [current_total_energy, current_avg_delay] = ...
    rismec_control2(N_slot, K, Nr, N_RIS, total_delta, N_blocks, qtz_RIS, weights, possible_angles, ...
        current_channel_ap_ris, current_channel_ap_ue, current_channel_ue_ris, proba);

    % Save current data
    total_energy(rr) = current_total_energy;
    avg_delay(rr) = current_avg_delay;

    simulation_time = toc;
    disp(['elapsed ', num2str(simulation_time), ' seconds.']);

end


