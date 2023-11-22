function [total_energy, avg_delay] = rismec_control2( ...
    N_slot, K, Nr, N_RIS, total_delta, N_blocks, qtz_RIS, weights, possible_angles, ...
    channel_ap_ris, channel_ap_ue, channel_ue_ris, proba)

    %% Test Probability Vector
    if sum(proba > 1)
         error('proba: Probability vector has more than one entry different than 1.');
    end


    %% Setup Parameters
    Msample = N_slot/2;                            % Slot to avoid for the final average
    perc = .9;                                     % Percentage of slot duration dedicated to offloading
    delta = total_delta .* perc;                   % Slot duration dedicated to offloading
    delta_signaling = (1 - perc) * total_delta;    % Slot duration dedicated to control signaling
    alpha = 0.5;                                   % Weight of user and network energy consumption (e.g. 0.5 --> holistic strategy)
    num_iter = 1;                                  % Number of iterations for possible sample average
    
    %% Communication Parameters
    Pt = 100e-3 * ones(K, 1);                      % Maximum user Tx power
    B = 100e6;                                     % Total bandwidth
    B_u = B / K;                                   % Total bandwidth UL per user
    F = 10^0.5 ;                                   % Noise figure
    N0 = 4e-21 * F;                                % Noise PSD
    p_AP_tx = 250e-3;                              % AP Tx power
    p_AP_on = 2.2;                                 % AP power ON
    
    % Power consumption of a RIS element, depending on the number of bits
    if qtz_RIS == 1
        p_bit = .5e-3;            
    elseif qtz_RIS == 2
        p_bit = 1e-3;
    elseif qtz_RIS == 3
        p_bit = 1.5e-3;
    end
    
    %% Computation Parameters 
    possible_f = (0:0.01:1) * 4.5e9;      % Possible CPU cycle frequencies
    min_poss_fm = 0.5e9;                  % Minimum needed CPU frequency
    kappa = 1e-27;                        % Effective switched capacitance of the processor
    x = 3 * ones(K, 1);                   % Conversion factor for the computation queue
    J = 10.^(-x);                         % Bits/CPU cycle
    
    %% Optimization Parameters and Requirements
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
    
    %% Simulation Results
    energy_tot = zeros(K, N_slot, length(V1)); energy_tot_mec = zeros(N_slot, length(V1));
    energy_tot_AP = energy_tot_mec;
    energy_tot_RIS = energy_tot_mec; 
    
    Z_tot = zeros(K, N_slot, length(V1));
    avg_delay_tot = Z_tot;
    tot_queues_ue = zeros(K, N_slot, length(V1)); 
    avg_delay = zeros(length(V1), num_iter);
    
    total_energy = zeros(length(V1), num_iter);
    total_energy_mec = total_energy;
    total_energy_AP = total_energy;
    total_energy_RIS = total_energy;

    %% Optimization Algorithm
    for iter = 1:num_iter

         % Prepare so save results in this iteration
         ue_comm_queue = zeros(K, N_slot, length(V1));      % local communication queue    
         es_remo_queue = zeros(K, N_slot, length(V1));      % remote computation queue
         es_comm_queue = zeros(K, N_slot, length(V1));      % virtual communication queue     
         
         % Loop over trade-off parameter V 
         for vv = 1:length(V1)
            V = V1(vv);       
            
            % Loop over the slots
            for nn = 1:N_slot
                
                % Shallow slot counting
                if mod(nn, 200) == 0
                    nn
                end

                % Toss a coin to evaluate possible errors for each UE
                toss = rand(K, 1);

                %% Phase 1: Signaling -- Initialization

                % Get current arrival arate
                A = A1(:, nn);     


                


                
                % Optimization function
                [power_up, rate_up, freq_MEH, power_mec, current_v, power_RIS, opt_beam] = ...
                optimization_RIS(K, B_u,...
                    ue_comm_queue(:, nn, vv), es_remo_queue(:, nn, vv), es_comm_queue(:, nn, vv), qtz_RIS, N_RIS, N_blocks,...
                    channel_ue_ris(:, :, nn), channel_ap_ris(:, :, :, nn), channel_ap_ue(:, :, nn),...
                    N0, V, Pt, possible_f, J, delta, p_bit, kappa, alpha, possible_angles, weights, Nr);
                
                % Update physical queues          
                ue_comm_queue(:, nn+1, vv) = max(0, ue_comm_queue(:,nn,vv) - delta * rate_up) + A;
                es_remo_queue(:,nn+1,vv) = max(0, es_remo_queue(:,nn,vv) - freq_MEH * delta .* J) + min(ue_comm_queue(:, nn, vv), delta * rate_up);
                
                % Update virtual queues
                es_comm_queue(:, nn+1, vv) = max(0, es_comm_queue(:, nn, vv) + ue_comm_queue(:, nn+1, vv) + es_remo_queue(:, nn+1, vv) - Q_avg);
                Z_tot(:, nn, vv) = Z_tot(:, nn, vv) + es_comm_queue(:, nn, vv); 
                
                % Energy donsumption and delay updates 
                energy_tot(:, nn, vv) = power_up * delta + Pt .* delta_signaling;
                energy_tot_AP(nn, vv) = delta * (p_AP_on * max(power_up>0)) + delta_signaling * (p_AP_on + p_AP_tx);
                energy_tot_mec(nn, vv) = delta * power_mec + delta_signaling * kappa * min_poss_fm^3;
                energy_tot_RIS(nn, vv) = delta * power_RIS + delta_signaling * p_bit * N_RIS;        
                tot_queues_ue(:, nn, vv) = ue_comm_queue(:, nn, vv) + es_remo_queue(:, nn, vv);
                avg_delay_tot(:, nn, vv) = avg_delay_tot(:, nn, vv) + mean(tot_queues_ue(:, nn - min(nn, 1e3) + 1:nn, vv), 2) ./ A_avg * delta * 1e3;
            end  
            
            total_energy(vv, iter) = mean(sum(energy_tot(:, Msample:N_slot, vv), 1), 2);
            total_energy_AP(vv, iter) = mean(energy_tot_AP(Msample:N_slot, vv));
            total_energy_mec(vv, iter) = mean(energy_tot_mec(Msample:N_slot, vv));
            total_energy_RIS(vv, iter) = mean(energy_tot_RIS(Msample:N_slot, vv));
            avg_delay(vv, iter) = max(mean(tot_queues_ue(:, Msample:N_slot, vv), 2) ./ A_avg * total_delta * 1e3);
            
         end
    end
end

