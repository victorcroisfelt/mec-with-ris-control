function [energy_tot, avg_delay_tot] = rismec_control( ...
    N_slot, K, N, f_max, tau, perc, N_blocks, qtz_RIS, weights, possible_angles, ...
    channel_ap_ris, channel_ap_ue, channel_ue_ris, error_prob, ...
    tti_time, conf_codebook_size, N_pilot, f_ra, tau_ra)

    %% Test Probability Vector
    if (length(error_prob) ~= 4)
        error('probas: probability vector has not 4 entries. You should have an entry for each error type, given in the order: INI-U, INI-R, SET-U, SET-R');
    end
    if any(error_prob > 1)
         error('probas: probability vector has more than one entry different than 1.');
    end
   
    %% Setup Parameters
    Msample = N_slot/10;                            % Slot to avoid for the final average                             
    tau_pay = tau .* perc;                         % Slot duration dedicated to offloading
    tau_ctl = (1 - perc) * tau;                    % Slot duration dedicated to control signaling
    sigma = 0.5;                                   % Weight of user and network energy consumption (e.g. 0.5 --> holistic strategy)
    num_iter = 1;                                  % Number of iterations for possible sample average
    
    %% Communication Parameters
    Pt = 100e-3 * ones(K, 1);                      % Maximum user Tx power [W]
    B = 100e6;                                     % Total bandwidth [Hz]
    B_u = B / K;                                   % Total bandwidth UL per user [Hz]
    F = 10^0.5 ;                                   % Noise figure
    N0 = 4e-21 * F;                                % Noise PSD
    p_AP_tx = 250e-3;                              % AP Tx power [W]
    p_AP_on = 2.2;                                 % AP power ON [W]
    
    % Power consumption of a RIS element, depending on the number of bits
    if qtz_RIS == 1
        p_bit = .5e-3;            
    elseif qtz_RIS == 2
        p_bit = 1e-3;
    elseif qtz_RIS == 3
        p_bit = 1.5e-3;
    end
    
    %% Computation Parameters 
    possible_f = (0:0.01:1) * f_max;      % Possible CPU cycle frequencies [Hz]    
    kappa = 1e-27;                        % Effective switched capacitance of the processor
    x = 3 * ones(K, 1);                   % Conversion factor for the computation queue
    J = 10.^(-x);                         % Bits/CPU cycle
    
    %% Optimization Parameters and Requirements
    A_avg = 200;                  % Average arrival rate    
    A = poissrnd(A_avg, K, N_slot);  % N_slot + 1, because it needs to 

    D_avg = 300e-3 * ones(K, 1);                % Delay constraint
    Q_avg = D_avg .* A_avg / tau;      % Average queue length corresponding to desired average delay
    
    %%% the vector V1 can be represented also by one value, usually the last one = the one for which the trade-off is optimal 
    %%% the vector is necessary only for the Energy/Delay trade-off
    % V1 = [1e8; 1e9; 2e9; 3e10];                % Lyapunov trade-off parameter
    V1 = 4e11;
    
    %% Simulation Results

    % Energy
    energy_tot_ues = zeros(K, N_slot, length(V1)); 
    energy_tot_es = zeros(N_slot, length(V1));
    energy_tot_ap = energy_tot_es;
    energy_tot_ris = energy_tot_es; 
    energy_tot = zeros(N_slot, length(V1));
   
    total_avg_energy = zeros(length(V1), num_iter);

    es_virtual_queue_total = zeros(K, N_slot, length(V1));
    avg_delay_tot = zeros(K, N_slot, length(V1));
    tot_queues_ue = zeros(K, N_slot, length(V1)); 
    avg_delay = zeros(length(V1), num_iter);

    %% Optimization Algorithm
    for iter = 1:num_iter

         % Prepare to save results in this iteration
         ue_comm_queue = zeros(K, N_slot, length(V1));         % UE's local communication queue    
         es_comp_queue = zeros(K, N_slot, length(V1));         % ES' computation queue
         es_comm_queue = zeros(K, N_slot, length(V1));     % UE's communication queue known by the ES         
         es_virtual_queue = zeros(K, N_slot, length(V1));      % ES' virtual communication queue (estimated)     

         % Prepare to save previous values
         prev_ue_rate = zeros(K, 1);            % previous UE's rate
         prev_ue_power = zeros(K, 1);           % previous UE's power     
         prev_ris_config = zeros(N, 1);     % previous RIS' configuration
         prev_freq_MEH = zeros(K, 1);           % previous UE's power     

         % Loop over trade-off parameter V 
         for vv = 1:length(V1)

            % Get current trade-off parameter
            V = V1(vv);       

            % Toss several coins to evaluate the packets that were losts
            tosses_ini_u = rand(N_slot, K);
            tosses_ini_r = rand(N_slot, 1);
            tosses_set_u = rand(N_slot, K);
            tosses_set_r = rand(N_slot, 1);

            % Compute packet lost
            pkts_lost_ini_u = tosses_ini_u < error_prob(1);
            pkts_lost_ini_r = tosses_ini_r < error_prob(2);
            pkts_lost_set_u = tosses_set_u < error_prob(3);
            pkts_lost_set_r = tosses_set_r < error_prob(4);    

            % Pre-allocation of data in the comm queue during first slot            
            
            
            % Loop over the slots
            for nn = 1:N_slot
                
                % Shallow slot counting
                if mod(nn, 50) == 0
                    fprintf('\t\titeration %03d/%03d\n', nn, N_slot);
                end

                % Get current, true channels
                true_channel_ue_ris = channel_ue_ris(:, :, nn);
                true_channel_ap_ris = channel_ap_ris(:, :, nn);
                true_channel_ap_ue = channel_ap_ue(:, :, nn);

                % Get current losses
                pkt_lost_ini_u = pkts_lost_ini_u(nn, :);
                pkt_lost_ini_r = pkts_lost_ini_r(nn);
                pkt_lost_set_u = pkts_lost_set_u(nn, :);
                pkt_lost_set_r = pkts_lost_set_r(nn);
              
                %% Phase 1: Signaling -- Initialization

                if nn == 1 
                    % On the first slot only the coomunication needs to be updated                    
                    % Update UE local comm queue
                    ue_comm_queue(:, nn, vv) = A(:, nn); 
                    % Update the ES estimated comm queue if the INI-U is
                    % correctly received
                    es_comm_queue(~pkt_lost_ini_u, nn, vv) = A(~pkt_lost_ini_u, nn, vv);
                    
                    % Update virtual queues
                    es_virtual_queue(:, nn, vv) = max(0, es_comm_queue(:, nn, vv) + es_comp_queue(:, nn, vv) - Q_avg);                                    
                    es_virtual_queue_total(:, nn, vv) = es_virtual_queue(:, nn, vv);    %% useless
                    
                else
                    % Update UE local comm queue
                    ue_comm_queue(:, nn, vv) = max(0, ue_comm_queue(:, nn-1, vv) - tau_pay * prev_ue_rate) + A(:, nn);
                    % Update the ES estimated comm queue 
                    %%% Check if INI-U was lost
                    % if the INI-U is correctly received                    
                    es_comm_queue(~pkt_lost_ini_u, nn, vv) = ue_comm_queue(~pkt_lost_ini_u, nn, vv);
                    % if the INI-U is not decoded
                    es_comm_queue(pkt_lost_ini_u, nn, vv) = max(0, ue_comm_queue(pkt_lost_ini_u, nn-1, vv) - tau_pay * prev_ue_rate(pkt_lost_ini_u));
                                                        
                    % Update ES (from its perpective)
                    es_comp_queue(:, nn, vv) = max(0, es_comp_queue(:, nn-1, vv) - prev_freq_MEH * tau_pay .* J) + min(es_comm_queue(:, nn-1, vv), tau_pay * prev_ue_rate);                    

                    % Update virtual queues
                    es_virtual_queue(:, nn, vv) = max(0, es_virtual_queue(:, nn-1, vv) + es_comm_queue(:, nn, vv) + es_comp_queue(:, nn, vv) - Q_avg);                                    
                    es_virtual_queue_total(:, nn, vv) = es_virtual_queue_total(:, nn, vv) + es_virtual_queue(:, nn, vv);    %% useless
                end 

                %%% Check if INI-R was lost 
                % moved after algorithmic phase to keep the ES freq optimized
                
                
                %% Phase 2: Algorithmic

                % Channel estimation
                hat_channel_ue_ris = true_channel_ue_ris;
                hat_channel_ap_ris = true_channel_ap_ris;
                hat_channel_ap_ue = true_channel_ap_ue;

                % Optimization 
                [ue_power, ue_rate, freq_MEH, power_mec, ris_config, power_RIS, ap_opt_beam] = ...
                optimization_RIS(K, B_u, ...
                    es_comm_queue(:, nn, vv), es_comp_queue(:, nn, vv), es_virtual_queue(:, nn, vv), ...
                    qtz_RIS, N, N_blocks, ...
                    hat_channel_ue_ris, hat_channel_ap_ris, hat_channel_ap_ue,...
                    N0, V, Pt, possible_f, J, tau_pay, p_bit, kappa, sigma, possible_angles, weights);


                %% Phase 3: Signaling -- Setup                
                %%% INI-R
                % If INI-R is lost no channel estimation can occur
                % we set the rates and powers to zero for this slot                
                if pkt_lost_ini_r
                    ue_rate = zeros(K, 1);
                    ue_power = zeros(K, 1);
                end                   

                %%% SET-U
                % If SET-U is lost no power and rate info is delivered to the UE
                % we override with the previous rates and values   
                ue_rate(pkt_lost_set_u) = prev_ue_rate(pkt_lost_set_u);
                ue_power(pkt_lost_set_u) = prev_ue_power(pkt_lost_set_u);           

                % Capacity computation 
                % (equivalent to evaluating the rate with true channels)
                %%% SET-R
                % If SET-R is lost the RISC loads the previous configuration                
                if pkt_lost_set_r
                    ue_capacity = eval_rate(true_channel_ue_ris, true_channel_ap_ris, true_channel_ap_ue, N0, B_u, prev_ris_config, weights(:, ap_opt_beam), ue_power);
                else
                    ue_capacity = eval_rate(true_channel_ue_ris, true_channel_ap_ris, true_channel_ap_ue, N0, B_u, ris_config, weights(:, ap_opt_beam), ue_power);
                end

                %% Phase 4: Payload                  
                %%% Test rate versus capacity
                ue_rate(ue_rate > ue_capacity) = 0;                    

                %% Update metrics                
                % Update energy consumption
                energy_tot_ues(:, nn, vv) = tau_pay * ue_power + Pt .* (tti_time + tti_time * conf_codebook_size * N_pilot);
                energy_tot_ap(nn, vv) = 2 * tti_time * (p_AP_on + p_AP_tx) * (K + 1);
                energy_tot_es(nn, vv) = tau_pay * power_mec + tau_ra * kappa * f_ra^3;
                energy_tot_ris(nn, vv) = tau_pay * power_RIS + (tau_ctl - tau_ra) * p_bit * N;  

                energy_tot(nn, vv) = sigma * sum(energy_tot_ues(:, nn, vv)) + (1-sigma) * (energy_tot_ris(nn, vv) + energy_tot_ap(nn, vv) + energy_tot_es(nn, vv));
                
                % Update delay
                tot_queues_ue(:, nn, vv) = ue_comm_queue(:, nn, vv) + es_comp_queue(:, nn, vv);
                avg_delay_tot(:, nn, vv) = avg_delay_tot(:, nn, vv) + mean(tot_queues_ue(:, nn - min(nn, 1e3) + 1:nn, vv), 2) ./ A_avg * tau_pay * 1e3;
            
                % Override previous values
                prev_ue_rate = ue_rate;
                prev_ue_power = ue_power;
                prev_ris_config = ris_config;
                prev_freq_MEH = freq_MEH;

            end  

           
%             
%             total_avg_energy(vv, iter) = mean(energy_tot(Msample:N_slot, vv), 1);
%             avg_delay(vv, iter) = max(mean(tot_queues_ue(:, Msample:N_slot, vv), 2) ./ A_avg * tau * 1e3);
            
         end

          energy_tot = squeeze(min(energy_tot, [], 2));
          avg_delay_tot = squeeze(min(avg_delay_tot, [], 3));
    end
end

