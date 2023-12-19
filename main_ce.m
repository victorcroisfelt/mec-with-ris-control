clear all; clc; rng('default');

%% Load Channel and Simulation Parameters 
file_name = 'data/channel/N64_M8_K4.mat';
load(file_name)

disp('----- Main Parameters -----')
disp(['N_slot = ', num2str(N_slot)])
disp(['K = ', num2str(K)])
disp(['N = ', num2str(N)])
disp(['M = ', num2str(M)])
disp(['num_setups = ', num2str(num_setups)])
disp('---------------------------')


% Load receive combining (CHECK THE CORRESPONDING N_AP)
load(['data/combining/w_candidates_' num2str(M) '_antennas.mat'], 'weights')

% Number of bits used to quantize the phase of the RIS
qtz_RIS = 1; 

% Max frequency available at the ES
f_max = 4.5e9;      % [Hz]

f_ra = 0.5e9;


% Define the RIS switching time
tau_ris = 0; % 1/14 * 10e-3

% Define kind of RIS-CC
ris_cc = 'ib-cc';
% ris_cc = 'ob-cc';

%% discount
discount_vec = 0.9; 0.1:0.2:1;

for discount = discount_vec

%% For N_blocks
% Number of blocks to reduce complexity (e.g. if N_RIS = 100 and N_blocks =
% 100/4 = 25, 4 elements together are optimized)
N_blocks_vec = [N/2, N];     
% be very careful in choosing the right divider respect to N_RIS

for N_blocks = N_blocks_vec

%% CE

% Number of pilots per TTI
N_p = 1;

% CE
conf_codebook_size = N + 0:16:128;

% Get length of the probability vector
num_points = length(conf_codebook_size);


%% Computing time
% Define range of total slot duration
tau = 100*1e-3;

% Define TTI time
constant = 1;
tti_time = constant * 1/14 * 10e-3;

[tau_overhead_range, perc_range, tau_sig, tau_ce, tau_ra] = compute_overhead(tti_time, tau_ris, tau, N_p, N_blocks, conf_codebook_size, qtz_RIS, M, f_ra, K, ris_cc);


%% Arrival rate and Lyapunov
avg_arrival_rate = [1, 100, 200, 500]*1e3;    % bit/s arrival rate corrisponding to 500 and 10000 average arrival bit at 100 ms of frame
num_arrival = length(avg_arrival_rate);

lyapunov_tradeoff = [1e8];   % Lyapunov trade-off parameters
num_lyapunov = length(lyapunov_tradeoff);

%% Prepare to save simulation results
total_energy = zeros(num_arrival, num_setups, num_points, N_slot, num_lyapunov);
avg_delay = zeros(num_arrival, num_setups, num_points, K, N_slot, num_lyapunov);

%% Simulation
disp('------- Simulation --------')

% Go through different ris-cc
%%% To be implemented

% Go through the possible errors
% In the first case no error, in the second the error is loaded considering
% the performance of the users according to the paper
for error_type = 0

   
% Go through all setups
for rr = 1:num_setups

    fprintf('setup = %02d/%02d\n', rr, num_setups)
    tic;

    % Get current channels
    current_channel_ris_ap = squeeze(overall_channel_ris_ap(rr, :, :, :, :));
    current_channel_ue_ap = squeeze(overall_channel_ue_ap(rr, :, :, :));
    current_channel_ue_ris = squeeze(overall_channel_ue_ris(rr, :, :, :));

    %%%  Error probability computation
    % possibility of computing the error probability depending on the setup  
    if error_type == 0
        error_prob = zeros(4,1);
    else
        %%%% not implemented yet
        error_prob = zeros(4,1);
    end

    %%% Go through different arrival time
    for aa = 1:num_arrival
        avg_arrival = avg_arrival_rate(aa);
        fprintf('\tarrival %02d/%02d\n', aa, num_arrival)

        % Go through all taus
        for nn = 1:num_points
           fprintf('\t\tpoint = %02d/%02d\n', nn, num_points)   
            
            % Recompute percentage of slot duration dedicated to offloading
            perc = perc_range(nn);

            % recompute ce codebook
            conf_codebook_size_now = conf_codebook_size(nn);
    
            % Optmize control
            [total_energy(aa, rr, nn, :, :), avg_delay(aa, rr, nn, :, :, :)] = ...
            rismec_control_ce(f_max, tau, perc, ...
                N_blocks, qtz_RIS, weights, ...
                current_channel_ris_ap, current_channel_ue_ap, current_channel_ue_ris, ...
                error_prob, tti_time, conf_codebook_size_now, N_p, f_ra, tau_ra, avg_arrival, lyapunov_tradeoff, ...
                discount, true);
 
        end
    end
    
    simulation_time = toc;
    disp(['elapsed ', num2str(simulation_time), ' seconds.']);
end

% Save results
file_name = ['data/results/ce/N', num2str(N), '_M', num2str(M), '_K', num2str(K), '_N_block', num2str(N_blocks), '_discount', num2str(discount, '%1.1f'), '.mat'];
save(file_name)

end

end
end
