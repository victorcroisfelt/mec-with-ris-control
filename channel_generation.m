clear all; clc; rand('state',0); randn('state',0) 

%% Parameters 
N_slot = 3e3;   % Number of slots                                  
N_RIS = 64;     % Number of RIS elements !!! NEEDS TO BE A PERFECT SQUARE FOR CHANNELS GENERATION
Nt = 1;         % Number of TX antennas (user) 
Nr = 32;        % Number of RX antennas, possibilities: 8, 16, 32 (AP)
alpha_dir = 3;  % FSPL exponent of the direct link
f = 2e9;        % Frequency

%% Simulation Parameters
K = 2;                  % Number of users (1, 2 or 4) 
dist_ap_ris = 100;      % Distance AP to RIS in meters
angl_ap_ris = 45;       % Angle AP to RIS in degrees
dist_minimum = 10;      % Minimum distance between RIS and UE
dist_maximum = 100;     % Maximum distance between RIS and UE
num_realizations = 10;  % Number of realizations

%% Simulation

% Prepare to save simulation results
overall_channel_ue_ris = zeros(num_realizations, N_RIS, K, N_slot); 
overall_channel_ap_ris = zeros(num_realizations, Nr, N_RIS, K, N_slot); 
overall_channel_ap_ue = zeros(num_realizations, Nr, K, N_slot);

% Go through all realizations
for rr = 1:num_realizations
    
    disp(['realization = ', num2str(rr)])

    % Go through each user
    for qq = 1:K
    
        % Generate channels
        [Hdirt, H1t, H2t] = channel_mat_RIS( ...
            Nt, Nr, N_RIS, ...
            dist_ap_ris, angl_ap_ris, ...
            dist_minimum, dist_maximum, ...
            f, N_slot, alpha_dir ...
            );
        
        % Get channels
        overall_channel_ue_ris(rr, :, qq, :) = cell2mat(H1t);
        overall_channel_ap_ris(rr, :, :, qq, :) = reshape(cell2mat(H2t), [Nr, N_RIS, N_slot]);
        overall_channel_ap_ue(rr, :, qq, :) = cell2mat(Hdirt);
    end

end

% Save generated channels
file_name = ['data/channel/K', num2str(K), '.mat'];
save(file_name)
