clear all; clc; rand('state',0); randn('state',0) 

%% Parameters 
N_slot = 3e3;   % Number of slots                                  
N_RIS = 64;     % Number of RIS elements !!! NEEDS TO BE A PERFECT SQUARE FOR CHANNELS GENERATION
Nt = 1;         % Number of TX antennas
Nr = 32;        % Number of RX antennas, possibilities: 8, 16, 32

% due to the particular geometric construction of channels, be carefull to
% change the number of UEs without changing the positions in terms of distances!
K = 4;                                    % Number of Users (1, 2 or 4) 
Dvec = [500, 400, 600, 300];              % TX-RX distance
dist_ris_vec = [40, 30, 50, 40];          % RIS distance from TX
lt_vec = [50, 40, 50, 30];                % TX position
lr = 50;                                  % RX position

rice_fac = 3;           % Rician factor
alpha_dir = 3;          % FSPL exponent of the direct link
f = 2e9;                % Frequency
overall_channel_user_RIS = zeros(N_RIS, K, N_slot); 
overall_channel_RIS_AP = zeros(Nr, N_RIS, K, N_slot); 
overall_channel_AP = zeros(Nr, K, N_slot);
for qq = 1:K
    D = Dvec(qq); dist_ris = dist_ris_vec(qq); lt = lt_vec(qq);
    [Hdirt, H1t, H2t] = channel_mat_RIS(Nt, Nr, N_RIS, lt, lr, D, N_slot, rice_fac, f, dist_ris, alpha_dir);
    overall_channel_user_RIS(:, qq, :) = cell2mat(H1t);
    overall_channel_RIS_AP(:, :, qq, :) = reshape(cell2mat(H2t), [Nr, N_RIS, N_slot]);
    overall_channel_AP(:, qq, :) = cell2mat(Hdirt);
end