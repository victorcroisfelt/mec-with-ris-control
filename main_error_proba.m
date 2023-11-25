clear all; clc; rng('default'); rng(0);

%% Load Channel and Simulation Parameters 
file_name = 'data/channel/N64_M8_K2.mat';
load(file_name)

disp('----- Main Parameters -----')
disp(['N_slot = ', num2str(N_slot)])
disp(['K = ', num2str(K)])
disp(['N = ', num2str(N)])
disp(['M = ', num2str(M)])
disp(['num_setups = ', num2str(num_setups)])
disp('---------------------------')

% Number of blocks to reduce complexity (e.g. if N_RIS = 100 and N_blocks =
% 100/4 = 25, 4 elements together are optimized)
N_blocks = N/2;     
% be very careful in choosing the right divider respect to N_RIS

% Load receive combining (CHECK THE CORRESPONDING N_AP)
load(['data/combining/w_candidates_' num2str(M) '_antennas.mat'], 'weights')
possible_angles = size(weights, 2);

% Number of bits used to quantize the phase of the RIS
qtz_RIS = 1; 

% Max frequency available at the ES
f_max = 4.5e9;      % [Hz]

f_ra = 0.5e9;

% Number of pilots per TTI
N_p = 1;

% Define the RIS switching time
tau_ris = 0; % 1/14 * 10e-3

% Define kind of RIS-CC
% ris_cc = 'ib-cc';
ris_cc = 'ob-cc';

%% Computing Time
% Define total slot duration
tau = 100e-3;

% Define TTI time
constant = 1;
tti_time = constant * 1/14 * 10e-3;

conf_codebook_size = N_blocks + 1;

[tau_overhead, perc, tau_sig, tau_ce, tau_ra] = compute_overhead(tti_time, tau_ris, tau, N_p, N_blocks, conf_codebook_size, qtz_RIS, M, f_ra, K, ris_cc);

%% Simulation Parameters
% Define error type
% 1: INI-U
% 2: INI-R
% 3: SET-U
% 4: SET-R

%error_type = 'ini-r';
%error_type = 'set-u';
%error_type = 'set-r';

% Probability vector to lose of the packets:
prob_vector = [0, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99];

% Get length of the probability vector
num_points = length(prob_vector);

%% Loop through error type

for error = 0:4
    switch error
        case 0
            error_type = 'ini-u';
            error_prob_matrix = [prob_vector.', zeros(num_points, 1), zeros(num_points, 1), zeros(num_points, 1)];
        case 1
            error_type = 'ini-r';
            error_prob_matrix = [zeros(num_points, 1), prob_vector.', zeros(num_points, 1),zeros(num_points, 1)];
        case 2
            error_type = 'set-u';
            error_prob_matrix = [zeros(num_points, 1), zeros(num_points, 1), prob_vector.', zeros(num_points, 1)];
        case 3
            error_type = 'set-r';
            error_prob_matrix =  [zeros(num_points, 1), zeros(num_points, 1), zeros(num_points, 1), prob_vector.'];     
   end

        
    disp(['----------', error_type, '------------'])

%% Prepare to save simulation results
total_energy = zeros(num_setups, num_points, N_slot);
avg_delay = zeros(num_setups, num_points, K, N_slot);

%% Simulation
disp('------- Simulation --------')

% Go through all setups
for rr = 1:num_setups

    fprintf('setup = %02d/%02d\n', rr, num_setups)
    tic;

    % Get current channels
    current_channel_ris_ap = squeeze(overall_channel_ris_ap(rr, :, :, :, :));
    current_channel_ue_ap = squeeze(overall_channel_ue_ap(rr, :, :, :));
    current_channel_ue_ris = squeeze(overall_channel_ue_ris(rr, :, :, :));

    % Go through all probability points
    parfor nn = 1:num_points       
        fprintf('\tpoint = %02d/%02d\n', nn, num_points)

        
        
        % Optmize control
        [total_energy(rr, nn, :), avg_delay(rr, nn, :, :)] = ...
        rismec_control(N_slot, K, N, ...
            f_max, tau, perc, ...
            N_blocks, qtz_RIS, weights, possible_angles, ...
            current_channel_ris_ap, current_channel_ue_ap, current_channel_ue_ris, ...
            error_prob_matrix(nn, :), tti_time, conf_codebook_size, N_p, f_ra, tau_ra);

    end

    simulation_time = toc;
    disp(['elapsed ', num2str(simulation_time), ' seconds.']);

end

% Save results
file_name = ['data/results/', error_type, '/N', num2str(N), '_M', num2str(M), '_K', num2str(K), '.mat'];
save(file_name)

end