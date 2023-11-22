function [total_energy, avg_delay] = rismec_control( ...
    N_slot, K, Nr, N_RIS, total_delta, perc, N_blocks, qtz_RIS, weights, possible_angles, ...
    channel_ap_ris, channel_ap_ue, channel_ue_ris, proba, error_type)

%     %% Test Probability Vector
%     if (length(probas) ~= 4)
%         error('probas: probability vector has not 4 entries. You should have an entry for each error type, given in the order: INI-U, INI-R, SET-U, SET-R');
%     end
% 
%     if sum(probas > 1)
%          error('probas: probability vector has more than one entry different than 1.');
%     end
% 
%     if sum(probas == 0) == 4
%         error('probas: all entries of the probability vector are 0.');
%     end
% 
%     % Get error index 
%     error_index = find(probas);
% 
%     % Define error type
%     if error_index == 1
%         error_type = 'ini-u';
%     elseif error_index == 2
%         error_type = 'ini-r';
%     elseif error_index == 3
%         error_type = 'set-u';
%     else
%         error_type = 'set-r';
%     end
   
    %% Setup Parameters
    Msample = N_slot/2;                            % Slot to avoid for the final average                             
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

    % Energy
    energy_tot = zeros(K, N_slot, length(V1)); 
    energy_tot_mec = zeros(N_slot, length(V1));
    energy_tot_AP = energy_tot_mec;
    energy_tot_RIS = energy_tot_mec; 
   
    total_energy = zeros(length(V1), num_iter);
    total_energy_mec = total_energy;
    total_energy_AP = total_energy;
    total_energy_RIS = total_energy;

    es_comm_queue_total = zeros(K, N_slot, length(V1));
    avg_delay_tot = zeros(K, N_slot, length(V1));
    tot_queues_ue = zeros(K, N_slot, length(V1)); 
    avg_delay = zeros(length(V1), num_iter);

    %% Optimization Algorithm
    for iter = 1:num_iter

         % Prepare to save results in this iteration
         ue_comm_queue = zeros(K, N_slot, length(V1));      % UE's local communication queue    
         es_remo_queue = zeros(K, N_slot, length(V1));      % ES' remote computation queue
         es_comm_queue = zeros(K, N_slot, length(V1));      % ES' virtual communication queue     
         
         % Prepare to save previous values
         prev_ue_rate = zeros(K, 1);            % previous UE's rate
         prev_ue_power = zeros(K, 1);           % previous UE's power     
         prev_ris_config = zeros(N_RIS, 1);     % previous RIS' configuration

         % Loop over trade-off parameter V 
         for vv = 1:length(V1)

            % Get current trade-off parameter
            V = V1(vv);       

            % Toss several coins to evaluate the packets that were losts
            if strcmp(error_type, 'ini-u') || strcmp(error_type, 'set-u')
                tosses = rand(N_slot, K);
            else 
                tosses = rand(N_slot, 1);
            end
            
            % Loop over the slots
            for nn = 1:N_slot
                
                % Shallow slot counting
                if mod(nn, 200) == 0
                    nn
                end

                % Get current, true channels
                true_channel_ue_ris = channel_ue_ris(:, :, nn);
                true_channel_ap_ris = channel_ap_ris(:, :, :, nn);
                true_channel_ap_ue = channel_ap_ue(:, :, nn);

                % Get current toss
                if strcmp(error_type, 'ini-u') || strcmp(error_type, 'set-u')
                    toss = tosses(nn, :);
                else 
                    toss = tosses(nn);
                end

                %% Phase 1: Signaling -- Initialization

                % Get current arrival rate
                A = A1(:, nn);     

                % Check if INI-U was lost
                if strcmp(error_type, 'ini-u')
                    packet_loss = toss < proba;
                end

                % Check if INI-R was lost
                if strcmp(error_type, 'ini-r')
                    packet_loss = toss < proba;
                end
                
                %% Phase 2: Algorithmic

                % Channel estimation
                hat_channel_ue_ris = true_channel_ue_ris;
                hat_channel_ap_ris = true_channel_ap_ris;
                hat_channel_ap_ue = true_channel_ap_ue;

                % Optimization 
                [ue_power, ue_rate, freq_MEH, power_mec, ris_config, power_RIS, ap_beamforming] = ...
                optimization_RIS(K, B_u, ...
                ue_comm_queue(:, nn, vv), es_remo_queue(:, nn, vv), es_comm_queue(:, nn, vv), ...
                qtz_RIS, N_RIS, N_blocks, ...
                hat_channel_ue_ris, hat_channel_ap_ris, hat_channel_ap_ue,...
                N0, V, Pt, possible_f, J, delta, p_bit, kappa, alpha, possible_angles, weights, Nr);


                %% Phase 3: Signaling -- Setup

                % Check if SET-U was lost
                if strcmp(error_type, 'set-u')
                    packet_loss = toss < proba;
                end

                % Check if SET-R was lost
                if strcmp(error_type, 'set-r')
                    packet_loss = toss < proba;
                end

                % Capacity computation
                %%% SET-R
                % If SET-R is lost we use the previous configuration
                if strcmp(error_type, 'set-r')
                    [~, ue_capacity, ~, ~, ~, ~, ~] = ...
                    optimization_RIS(K, B_u, ...
                    ue_comm_queue(:, nn, vv), es_remo_queue(:, nn, vv), es_comm_queue(:, nn, vv), ...
                    qtz_RIS, N_RIS, N_blocks, ...
                    true_channel_ue_ris, true_channel_ap_ris, true_channel_ap_ue,...
                    N0, V, Pt, possible_f, J, delta, p_bit, kappa, alpha, possible_angles, weights, Nr, prev_ris_config);
                else
                    [~, ue_capacity, ~, ~, ~, ~, ~] = ...
                    optimization_RIS(K, B_u, ...
                    ue_comm_queue(:, nn, vv), es_remo_queue(:, nn, vv), es_comm_queue(:, nn, vv), ...
                    qtz_RIS, N_RIS, N_blocks, ...
                    true_channel_ue_ris, true_channel_ap_ris, true_channel_ap_ue,...
                    N0, V, Pt, possible_f, J, delta, p_bit, kappa, alpha, possible_angles, weights, Nr);
                end

                %% Phase 4: Payload
                
                % Do the optimization 
                %%% INI-R
                % If INI-R is lost we use the control configuration
                % we set the rates and powers to zero
                if strcmp(error_type, 'ini-r')
                    if packet_loss
                        ue_rate = zeros(K, 1);
                        ue_power = zeros(K, 1);
                    end   
                end                 

                %%% SET-U
                % If SET-U is lost we override with the previous rates
                % and values
                if strcmp(error_type, 'set-u')
                    for kk = 1:K
                        if packet_loss(kk)
                            ue_rate(kk) = prev_ue_rate(kk);
                            ue_power(kk) = prev_ue_power(kk);
                        end   
                    end
                end

                %%% Test rate versus capacity
                for kk = 1:K
                    if ue_rate(kk) > ue_capacity(kk)
                        ue_rate(kk) = 0;
                    end
                end
     
                % Update UE's queue
                %%% INI-U
                % If INI-U is lost we do not consider the knowledge of the
                % arrival rate when updating
                if strcmp(error_type, 'ini-u')
                    for kk = 1:K
                        if packet_loss(kk)
                            ue_comm_queue(kk, nn+1, vv) = max(0, ue_comm_queue(kk, nn, vv) - delta * ue_rate(kk));
                        else
                            ue_comm_queue(kk, nn+1, vv) = max(0, ue_comm_queue(kk, nn, vv) - delta * ue_rate(kk)) + A(kk);
                        end
                    end
                else
                    ue_comm_queue(:, nn+1, vv) = max(0, ue_comm_queue(:, nn, vv) - delta * ue_rate) + A;
                end

                % Update ES queues         
                es_remo_queue(:, nn+1, vv) = max(0, es_remo_queue(:, nn, vv) - freq_MEH * delta .* J) + min(ue_comm_queue(:, nn, vv), delta * ue_rate);
                es_comm_queue(:, nn+1, vv) = max(0, es_comm_queue(:, nn, vv) + ue_comm_queue(:, nn+1, vv) + es_remo_queue(:, nn+1, vv) - Q_avg);
                es_comm_queue_total(:, nn, vv) = es_comm_queue_total(:, nn, vv) + es_comm_queue(:, nn, vv); 
                
                % Update energy consumption
                energy_tot(:, nn, vv) = ue_power * delta + Pt .* delta_signaling;
                energy_tot_AP(nn, vv) = delta * (p_AP_on * max(ue_power>0)) + delta_signaling * (p_AP_on + p_AP_tx);
                energy_tot_mec(nn, vv) = delta * power_mec + delta_signaling * kappa * min_poss_fm^3;
                energy_tot_RIS(nn, vv) = delta * power_RIS + delta_signaling * p_bit * N_RIS;        
                
                % Update delay
                tot_queues_ue(:, nn, vv) = ue_comm_queue(:, nn, vv) + es_remo_queue(:, nn, vv);
                avg_delay_tot(:, nn, vv) = avg_delay_tot(:, nn, vv) + mean(tot_queues_ue(:, nn - min(nn, 1e3) + 1:nn, vv), 2) ./ A_avg * delta * 1e3;
            
                % Override previous values
                prev_ue_rate = ue_rate;
                prev_ue_power = ue_power;
                prev_ris_config = ris_config;

            end  
            
            total_energy(vv, iter) = mean(sum(energy_tot(:, Msample:N_slot, vv), 1), 2);
            total_energy_AP(vv, iter) = mean(energy_tot_AP(Msample:N_slot, vv));
            total_energy_mec(vv, iter) = mean(energy_tot_mec(Msample:N_slot, vv));
            total_energy_RIS(vv, iter) = mean(energy_tot_RIS(Msample:N_slot, vv));
            avg_delay(vv, iter) = max(mean(tot_queues_ue(:, Msample:N_slot, vv), 2) ./ A_avg * total_delta * 1e3);
            
         end
    end
end

