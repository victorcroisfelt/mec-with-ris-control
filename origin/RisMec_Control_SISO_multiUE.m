
%% THIS SCRIPT IS VALID FOR SISO AND MULTIPLE UEs

clear all; clc; rand('state',0); randn('state',0)      

% Optimization parameters and data
N_slot = 3e3;                                  % Number of slots
Msample = N_slot/2;                            % slot to avoid for the final average
total_delta = 10e-3;                           % Total slot duration
perc = .9;                                     % Percentage of slot duration dedicated to offloading
delta = total_delta .* perc;                   % Slot duration dedicated to offloading
delta_signaling = (1 - perc) * total_delta;    % Slot duration dedicated to control signaling
alpha = 0.5;                                   % Weight of user and network energy consumption (e.g. 0.5 --> holistic strategy)
num_iter = 1;                                  % number of iterations for possible sample average

%% RIS parameters 
% Number of RIS elements !!! NEEDS TO BE A PERFECT SQUARE FOR CHANNELS GENERATION
N_RIS = 100;                   % e.g. tried also for 64, 81, 144, 
% Number of blocks to reduce complexity (e.g. if N_RIS = 100 and N_blocks = 100/4 = 25, 4 elements together are optimized)
N_blocks = N_RIS/1;            % be very careful in choosing the right divider respect to N_RIS
                               
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
Nr = 1;                 % Number of RX antennas 
rice_fac = 3;           % Rician factor
alpha_dir = 3;          % FSPL exponent of the direct link
f = 2e9;                % Frequency

% due to the particular geometric constraction of channels, be carefull to
% change the number of UEs without changing the positions in terms of distances!
K = 4;                                    % Number of Users(1, 2 or 4)
Dvec = [500, 400, 600, 300];              % TX-RX distance
dist_ris_vec = [40, 30, 50, 40];          % RIS distance from TX
lt_vec = [50, 40, 50, 30];                % TX position
lr = 50;                                  % RX position

overall_channel_user_RIS = zeros(N_RIS, K, N_slot); overall_channel_RIS_AP = zeros(K, N_RIS, N_slot); overall_channel_AP = zeros(K, N_slot);
for qq = 1:K
    D = Dvec(qq); dist_ris = dist_ris_vec(qq); lt = lt_vec(qq);
    [Hdirt, H1t, H2t] = channel_mat_RIS(Nt, Nr, N_RIS, lt, lr, D, N_slot, rice_fac, f, dist_ris, alpha_dir);
    overall_channel_user_RIS(:, qq, :) = cell2mat(H1t);
    overall_channel_RIS_AP(qq, :, :) = reshape(cell2mat(H2t), N_RIS, N_slot);
    overall_channel_AP(qq, :) = cell2mat(Hdirt);
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
energy_tot = zeros(K, N_slot, length(V1)); energy_tot_mec = zeros(N_slot, length(V1)); energy_tot_AP = energy_tot_mec; energy_tot_RIS = energy_tot_mec; 
Z_tot = zeros(K, N_slot, length(V1)); avg_delay_tot = Z_tot; tot_queues_ue = zeros(K, N_slot, length(V1)); 
total_energy = zeros(length(V1), num_iter); avg_delay = total_energy; total_energy_mec = total_energy; total_energy_AP = total_energy; total_energy_RIS = total_energy;

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
            [power_up, rate_up, freq_MEH, power_mec, current_v, power_RIS] = optimization_RIS(K, B_u,...
                Q_local(:, nn, vv), Q_MEH(:, nn, vv), Z(:, nn, vv), qtz_RIS, N_RIS, N_blocks,...
                overall_channel_user_RIS(:, :, nn), overall_channel_RIS_AP(:, :, nn), overall_channel_AP(:, nn),...
                N0, V, Pt, possible_f, J, delta, p_bit, kappa, alpha);
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



%% List of Functions %%%%---------------------------

%%% CHANNEL %%%
function [Hdir, H1, H2] = channel_mat_RIS(Nt, Nr, Nris, lt, lr, D, no_mat, Nslot, f, dist_ris, varargin)
% data
lambda = 3e8/f;     % Wavelength
dt = lambda/2;      % TX antenna space
dr = lambda/2;      % RX antenna space
dris = lambda/2;    % RIS element space
k = 2*pi/lambda;    % Wavenumber
% Geometrical placement, x, y and z axis 
% TX antenna array
tx_arr(1,:) = zeros(1,Nt); 
tx_arr(2,:) = (sort(0:Nt-1,'descend')-(Nt-1)/2)*dt+lt; 
tx_arr(3,:) = zeros(1,Nt); 
% RX antenna array
rx_arr(1,:) = D*ones(1,Nr);
rx_arr(2,:) = (sort(0:Nr-1,'descend')-(Nr-1)/2)*dr+lr; 
rx_arr(3,:) = zeros(1,Nr);
% RIS 
center = [dist_ris 0]; % RIS center position 
N1 = sqrt(Nris);
N2 = N1;                                    % Number of RIS elements in two dimensions N1 and N2
ris_pos = RISPosition(N1,N2,dris,center);   % RIS elements' coordinates 
a = repmat(ris_pos{1},N1,1);                % Placing RIS elements in proper coordinates
ris_arr(1,:) = a(:)';        
ris_arr(2,:) = zeros(1,Nris);
ris_arr(3,:) = repmat(ris_pos{2},1,N2); 

if isempty(varargin)                        % Load the FSPL of the direct link
    alpha = 2;
else
    alpha = varargin{1};
end

% direct TX-RX paths/channel matrix
for i1 = 1:Nr                                        % Distance between the TX and RX antennas                                           
    for j1 = 1:Nt
        d(i1,j1) = norm(rx_arr(:,i1)-tx_arr(:,j1));
    end 
end
Hdir_los = exp(-1i*k*d);                             % Direct link, LOS matrix exponents 
tx_rx_dist = sqrt(D^2+(lt-lr)^2);                    % TX-RX distance   
FSPL_dir = (lambda/(4*pi))^2/tx_rx_dist^alpha(1);    % Inversion of the FSPL of the direct link 
Hdir = Rician(Hdir_los,sqrt(FSPL_dir),no_mat,Nslot);     % Direct link channel matrix           

% indirect paths (TX-RIS-RX)
for l1 = 1:Nris                                      % Distance between the RIS elements and the RX antennas                                                            
    for r1 = 1:Nr  
        d2(r1,l1) = norm(rx_arr(:,r1)-ris_arr(:,l1)); 
    end
    for t1 = 1:Nt                                    % Distance between the RIS elements and the TX antennas                  
        d1(l1,t1) = norm(tx_arr(:,t1)-ris_arr(:,l1));   
    end
end

tx_ris_dist = sqrt(dist_ris^2+lt^2);                 % TX-RIS distance
ris_rx_dist = sqrt((D-dist_ris)^2+lr^2);             % RIS-RX distance   

FSPLindir = lambda^4/(256*pi^2)*...                  % Inversion of the FSPL of the indirect link 
           ((lt/tx_ris_dist+lr/ris_rx_dist)^2)*...
           1/(tx_ris_dist*ris_rx_dist)^2;

% TX-RIS channel matrix
H1_los = exp(-1i*k*d1);                             % TX-RIS link, LOS matrix exponents  
FSPL_1 = sqrt(FSPLindir);                           % FSPL of the indirect link is embedded in the TX-RIS channel matrix 
H1 = Rician(H1_los,FSPL_1,no_mat,Nslot);

% RIS-RX channel matrix
H2_los = exp(-1i*k*d2);                             % RIS-RX link, LOS matrix exponents 
FSPL_2 = 1;
H2 = Rician(H2_los,FSPL_2,no_mat,Nslot);
end

function pos = RISPosition(N1,N2,dist,center)                 % Determine positions of RIS elements
d1 = (0:N1-1)-(N1-1)/2;
d2 = (0:N2-1)-(N2-1)/2;
pos{1} = center(1)+d1*dist;
pos{2} = center(2)+d2*dist;
end 

function Hout = Rician(Hlos,FSPL,no_mat,K)                     % Create the Rician channel matices
Hlos = repmat(Hlos,no_mat,1);
Hnlos = sqrt(1/2)*(randn(size(Hlos))+1i*randn(size(Hlos)));
Htot = FSPL/sqrt(K+1)*(Hlos*sqrt(K)+Hnlos);
dim = size(Hlos,1)/no_mat;
for ind = 1:no_mat
   Hout{ind} = Htot((ind-1)*dim+1:ind*dim,:); 
end
end

%%% OPTIMIZATION %%%
function [optimal_power_up, optimal_rate_up, freq_MEH, power_mec, current_v, power_RIS] = optimization_RIS(K, B_u,...
    Q_local, Q_MEH, Z, qtz_RIS, N_RIS, N_blocks, channel_RIS, channel_RIS_AP, channel_AP,...
    N0, V, Pt, possible_f, J, delta, p_bit, kappa, alpha)

    % Rate Optimization
    Q_tilde_up = B_u * (Q_local - Q_MEH + Z);
    Q_tilde_comp = (Q_MEH + Z) .* J;
    m = [0:2^(qtz_RIS)-1]';
    possible_v = [0; exp(1i*(2*m*pi) / (2^(qtz_RIS)))];

    current_v = zeros(N_RIS, 1);
    current_v = Greedy_RIS(current_v, N_blocks, N_RIS, possible_v, K, channel_RIS, channel_RIS_AP,...
           channel_AP, N0, B_u, V, Q_tilde_up, p_bit, alpha);
   
    overall_channel_RIS_up = zeros(K, 1); channel_over_noise_up = zeros(K, 1);
    for kk = 1:K
        overall_channel_RIS_up(kk) = overall_channel_RIS_up(kk) + (channel_RIS_AP(kk, :) * diag(current_v) * channel_RIS(:, kk));
        channel_over_noise_up(kk, 1) = abs(channel_AP(kk) + overall_channel_RIS_up(kk)) .^2 / (N0*B_u);
        % Uplink 
        rate_max_up(kk) = min(B_u * log2(1 + channel_over_noise_up(kk) .* Pt(kk)), Q_local(kk) / delta);
        p_max(kk) = 1 ./ channel_over_noise_up(kk) .* (exp(rate_max_up(kk) * log(2) / B_u) - 1);
        optimal_power_up(kk,1) = min(max(0, Q_tilde_up(kk) / (alpha * V * log(2)) - 1 ./ (channel_over_noise_up(kk))), p_max(kk));          
    end
    optimal_power_up(Q_tilde_up<=0) = 0;
    optimal_rate_up = B_u * log2(1 + channel_over_noise_up .* optimal_power_up);
    power_RIS = p_bit * sum(abs(current_v));

    % CPU Scheduling Optimization
    objective_current = +Inf;
    for ff = 1:length(possible_f)
        Z1 = Z;
        f_tot = possible_f(ff);
        freq_MEH = zeros(K, 1);
        freq = f_tot;
        ii = 0;
        while freq>0
            ii = ii+1;
            [~, index] = max((Z1 + Q_MEH) .* J);
            index = index(1);
            freq_MEH(index) = min((Q_MEH(index)) ./ (delta * J(index)), freq);
            freq = freq - freq_MEH(index);
            Z1(index) = -Inf;
            if ii == K
                break
            end
        end
        freq_MEH(Q_tilde_comp<=0) = 0;
        objective = V * delta * (kappa * f_tot^3) * (1 - alpha) - sum(delta * freq_MEH .* Q_tilde_comp);
        if objective<objective_current
            objective_current = objective;
            freq_MEH_current = freq_MEH;
            f_opt = f_tot;
        else         
        end
    end
    freq_MEH = freq_MEH_current;
    power_mec = kappa * f_opt^3;
end

%%% GREEDY RIS OPTIMIZATION FUNCTION %%%
function current_v = Greedy_RIS(current_v, N_blocks, N_RIS, possible_v, K, channel_RIS, channel_RIS_AP,...
    channel_AP, N0, B_u, V, Q_tilde_up, p_bit, alpha)
    N_RIS_block = N_RIS / N_blocks;
    for ll = 1:N_blocks
        current_objective_greedy = +Inf;
        try_v = current_v;
        for vv = 1:length(possible_v)           
            overall_channel_RIS_up = zeros(K, 1); channel_over_noise_up = zeros(K, 1);
            try_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block) = possible_v(vv);
            for kk = 1:K
                overall_channel_RIS_up(kk) = overall_channel_RIS_up(kk) + (channel_RIS_AP(kk, :) * diag(try_v) * channel_RIS(:, kk));
                channel_over_noise_up(kk) = abs(channel_AP(kk) + overall_channel_RIS_up(kk)) .^2 / (N0 * B_u);
            end

            objective_greedy = V * p_bit * sum(abs(try_v)) * (1 - alpha) - sum(channel_over_noise_up .* max(Q_tilde_up, 0));
            
            if objective_greedy <= current_objective_greedy
                current_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block) = try_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block);
                current_objective_greedy = objective_greedy;        
            end
        end
    end     
end
