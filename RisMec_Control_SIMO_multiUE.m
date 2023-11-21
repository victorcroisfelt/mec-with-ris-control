
%% THIS SCRIPT IS VALID FOR SIMO AND MULTIPLE UEs
clear all; clc; rand('state',0); randn('state',0)      

% Optimization parameters and data
N_slot = 3e3;                                  % Number of slots
% it should be better to run for 3 or more for final resuts
Msample = N_slot/2;                            % Slot to avoid for the final average
total_delta = 10e-3;                           % Total slot duration
perc = .9;                                     % Percentage of slot duration dedicated to offloading
delta = total_delta .* perc;                   % Slot duration dedicated to offloading
delta_signaling = (1 - perc) * total_delta;    % Slot duration dedicated to control signaling
alpha = 0.5;                                   % Weight of user and network energy consumption (e.g. 0.5 --> holistic strategy)
num_iter = 1;                                  % Number of iterations for possible sample average

%% RIS parameters 
% Number of RIS elements !!! NEEDS TO BE A PERFECT SQUARE FOR CHANNELS GENERATION
N_RIS = 64;                   % e.g. tried also for 81, 100, 144, 
% Number of blocks to reduce complexity (e.g. if N_RIS = 100 and N_blocks = 100/4 = 25, 4 elements together are optimized)
N_blocks = N_RIS/2;            % be very careful in choosing the right divider respect to N_RIS

qtz_RIS = 1;                   % Number of bits used to quantize the phase
if qtz_RIS == 1
    p_bit = .5e-3;             % Power consumption of a RIS element, depending on the number of bits
elseif qtz_RIS == 2
    p_bit = 1e-3;
elseif qtz_RIS == 3
    p_bit = 1.5e-3;
end

%% Number of UEs and Channel generation (SISO)
Nt = 1;                 % Number of TX antennas
Nr = 32;                % Number of RX antennas, other possibilities: 8, 16
% be careful in choosing the right data repository:
% 'w_candidates_8_antennas.mat' or 'w_candidates_16_antennas.mat'
load('w_candidates_32_antennas.mat', 'weights')
possible_angles = size(weights,2);

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

%% Communication parameters
Pt = 100e-3 * ones(K, 1);                      % Maximum user Tx power
B = 100e6;                                     % Total bandwidth
B_u = B / K;                                   % Total bandwidth UL per user
F = 10^0.5 ;                                   % Noise figure
N0 = 4e-21 * F;                                % Noise PSD
p_AP_tx = 250e-3;                              % AP Tx power
p_AP_on = 2.2;                                 % AP power ON

%% Computation parameters 
possible_f = (0:0.01:1) * 4.5e9;      % Possible CPU cycle frequencies
min_poss_fm = 0.5e9;                  % Minimum needed CPU frequency
kappa = 1e-27;                        % Effective switched capacitance of the processor
x = 3 * ones(K, 1);                   % Conversion factor for the computation queue
J = 10.^(-x);                         % Bits/CPU cycle

%% Optimization parameters and requirements
A_avg = 200 * ones(K, 1);                  % Arrival rate
A1 = zeros(K, N_slot);
for kk = 1:K
    A1(kk,:) = poissrnd(A_avg(kk), 1, N_slot);
end
D_avg = 30e-3 * ones(K, 1);                % Delay constraint
Q_avg = D_avg .* A_avg / total_delta;      % Average queue length corresponding to desired average delay

%%% the vector V1 can be represented also by one value, usually the last one = the one for which the trade-off is optimal 
%%% the vector is necessary only for the Energy/Delay trade-off
V1 = [1e8; 1e9; 2e9; 3e10];                % Lyapunov trade-off parameter
V1 = 4e11;

%% Structures Inizialization
energy_tot = zeros(K, N_slot, length(V1)); energy_tot_mec = zeros(N_slot, length(V1));
energy_tot_AP = energy_tot_mec;
energy_tot_RIS = energy_tot_mec; 

Z_tot = zeros(K, N_slot, length(V1));
avg_delay_tot = Z_tot;
tot_queues_ue = zeros(K, N_slot, length(V1)); 
avg_delay = total_energy;

total_energy = zeros(length(V1), num_iter);
total_energy_mec = total_energy;
total_energy_AP = total_energy;
total_energy_RIS = total_energy;

%% Optimization Algorithm
for iter = 1:num_iter
     Q_local = zeros(K, N_slot, length(V1)); Q_MEH = Q_local;
     Z = zeros(K, N_slot, length(V1)); Y = Z; arrivals = Z;
     rand('state',iter); randn('state',iter)
     % Loop over trade-off parameter V 
     for vv = 1:length(V1)
        V = V1(vv);       
        % Loop over the Slots
        for nn = 1:N_slot
            % shallow slot counting
            if mod(nn, 200) == 0
                nn
            end
            A = A1(:, nn);            
            % Optimization Function
            [power_up, rate_up, freq_MEH, power_mec, current_v, power_RIS, opt_beam] = optimization_RIS(K, B_u,...
                Q_local(:, nn, vv), Q_MEH(:, nn, vv), Z(:, nn, vv), qtz_RIS, N_RIS, N_blocks,...
                overall_channel_user_RIS(:, :, nn), overall_channel_RIS_AP(:, :, :, nn), overall_channel_AP(:, :, nn),...
                N0, V, Pt, possible_f, J, delta, p_bit, kappa, alpha, possible_angles, weights, Nr);
            % Physical queues update          
            Q_local(:, nn+1, vv) = max(0, Q_local(:,nn,vv) - delta * rate_up) + A;
            Q_MEH(:,nn+1,vv) = max(0, Q_MEH(:,nn,vv) - freq_MEH * delta .* J) + min(Q_local(:, nn, vv), delta * rate_up);
            % Virtual queues update
            Z(:, nn+1, vv) = max(0, Z(:, nn, vv) + Q_local(:, nn+1, vv) + Q_MEH(:, nn+1, vv) - Q_avg);
            Z_tot(:, nn, vv) = Z_tot(:, nn, vv) + Z(:, nn, vv); 
            % Energy Consumption and delay Updates 
            energy_tot(:, nn, vv) = power_up * delta + Pt .* delta_signaling;
            energy_tot_AP(nn, vv) = delta * (p_AP_on * max(power_up>0)) + delta_signaling * (p_AP_on + p_AP_tx);
            energy_tot_mec(nn, vv) = delta * power_mec + delta_signaling * kappa * min_poss_fm^3;
            energy_tot_RIS(nn, vv) = delta * power_RIS + delta_signaling * p_bit * N_RIS;        
            tot_queues_ue(:, nn, vv) = Q_local(:, nn, vv) + Q_MEH(:, nn, vv);
            avg_delay_tot(:, nn, vv) = avg_delay_tot(:, nn, vv) + mean(tot_queues_ue(:, nn - min(nn, 1e3) + 1:nn, vv), 2) ./ A_avg * delta * 1e3;
        end  
        
        total_energy(vv, iter) = mean(sum(energy_tot(:, Msample:N_slot, vv), 1), 2);
        total_energy_AP(vv, iter) = mean(energy_tot_AP(Msample:N_slot, vv));
        total_energy_mec(vv, iter) = mean(energy_tot_mec(Msample:N_slot, vv));
        total_energy_RIS(vv, iter) = mean(energy_tot_RIS(Msample:N_slot, vv));
        avg_delay(vv, iter) = max(mean(tot_queues_ue(:, Msample:N_slot, vv), 2) ./ A_avg * total_delta * 1e3);
     end
end

%%% Plots to control the correct computation offloading procedure
%%% this plot is about the average delay which is the most important
%%% parameter to evaluate if the queues are converging = the delay constraint is satisfied
% = is under the threshold; futhermore, the delay has to increase for greater V1 values
figure
loglog(V1,avg_delay,'-s','Linewidth',4,'MarkerSize',12);
grid on
hold on
yline(1e3*D_avg,'-.','alpha',1,'Linewidth',3)
xlabel('Lyapunov trade-off parameter $V$','Interpreter','Latex')
ylabel('Average Delay','Interpreter','Latex')

%% ---------------------
%%% eventual average over multiple iterations
Z_tot = Z_tot/num_iter;
avg_delay_tot = avg_delay_tot/num_iter;

for vv = 1:length(V1)
    final_energy(vv) = mean(total_energy(vv,:));
    final_energy_AP(vv) = mean(total_energy_AP(vv,:));
    final_energy_mec(vv) = mean(total_energy_mec(vv,:));
    final_energy_RIS(vv) = mean(total_energy_RIS(vv,:));
    final_delay(vv) = mean(avg_delay(vv,:));
end

%% ----------------------
%%% control over the virtual queues
figure
for kk =1:K
    plot(Z_tot(kk,:,end),'-','Linewidth',4,'MarkerSize',12);
    hold on
end
grid on
ylabel('Virtual Queue','Interpreter','Latex')
xlabel('time slots','Interpreter','Latex')





