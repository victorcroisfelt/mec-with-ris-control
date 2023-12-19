function [optimal_power_ul, optimal_rate_up, freq_MEH, power_mec, opt_conf, power_RIS, opt_beam] = ...
    optimization_RIS_CE(K, B_u, Q_ue_comm, Q_es_proc, Q_virtual, ...
        qtz_RIS, N, N_blocks, direct_channel, reflected_channel, ...
        N0, V, Pt, possible_f, J, tau_pay, p_bit, kappa, alpha, bf_weights, varargin)

    % Rate Optimization parameters
    Q_tilde_up = B_u * (Q_ue_comm - Q_es_proc + Q_virtual);
    Q_tilde_comp = (Q_es_proc + Q_virtual) .* J;
    m = (0:2^(qtz_RIS)-1).';
    possible_confs = [0; exp(1i*(2*m*pi) / (2^(qtz_RIS)))];

    % Optimize the communication parameters
    current_conf = zeros(N, 1);
    [opt_conf, opt_beam] = greedy_optimization_ce(current_conf, N_blocks, N, possible_confs, ...
           direct_channel, reflected_channel, ...
           N0, B_u, V, Q_tilde_up, p_bit, alpha, bf_weights);
        
    if ~isempty(varargin)                        
        opt_conf = varargin{1};
    end

    % Compute optimal rate and power    
    channel_over_noise_ul = (abs(bf_weights(:, opt_beam)' * (direct_channel + squeeze(pagemtimes(permute(reflected_channel, [2, 3, 1]), opt_conf)))) .^2 / (N0 * B_u)).';       
    rate_max_up = min(B_u * log2(1 + channel_over_noise_ul .* Pt), Q_ue_comm / tau_pay);
    p_max = 1 ./ channel_over_noise_ul .* (exp(rate_max_up * log(2) / B_u) - 1);
    optimal_power_ul = min(max(0, Q_tilde_up / (alpha * V * log(2)) - 1 ./ (channel_over_noise_ul)), p_max);          
    % Put to zero power giving negative total normalized queue state
    optimal_power_ul(Q_tilde_up<=0) = 0;
    % Compute optimal rate
    optimal_rate_up = B_u * log2(1 + channel_over_noise_ul .* optimal_power_ul);
    power_RIS = p_bit * sum(abs(opt_conf));

    % CPU Scheduling Optimization
    objective_current = +Inf;
    for ff = 1:length(possible_f)
        Z1 = Q_virtual;
        f_tot = possible_f(ff);
        freq_MEH = zeros(K, 1);
        freq = f_tot;
        ii = 0;
        while freq>0
            ii = ii+1;
            [~, index] = max((Z1 + Q_es_proc) .* J);
            index = index(1);
            freq_MEH(index) = min(Q_es_proc(index) ./ (tau_pay * J(index)), freq);
            freq = freq - freq_MEH(index);
            Z1(index) = -Inf;
            if ii == K
                break
            end
        end
        freq_MEH(Q_tilde_comp<=0) = 0;
        objective = V * tau_pay * (kappa * f_tot^3) * (1 - alpha) - sum(tau_pay * freq_MEH .* Q_tilde_comp);
        if objective<objective_current
            objective_current = objective;
            freq_MEH_current = freq_MEH;
            f_opt = f_tot;        
        end
    end
    freq_MEH = freq_MEH_current;
    power_mec = kappa * f_opt^3;
end


%%% GREEDY OPTIMIZATION FUNCTION for RIS and COMBINER %%%
function [current_conf, opt_beam] = greedy_optimization_ce( ...
    current_conf, N_blocks, N, possible_confs, ...
    direct_channel, reflected_channel, ...
    N0, B_u, V, Q_tilde_up, p_bit, alpha, bf_weights)

    N_g = N / N_blocks;
    possible_angles = size(bf_weights, 2);

    % Checking for each block of RIS elements
    for ll = 1:N_blocks
        current_objective_greedy = +Inf; current_obj_prov = +Inf;
        try_conf = current_conf;
        for vv = 1:length(possible_confs)           
            try_conf((ll - 1) * N_g + 1 : ll * N_g) = possible_confs(vv);
            
            %%% looking for the optimal AP combiner
            for bb = 1:possible_angles
                % Compute C/N with combiner and configuration under test
                % transposition to obtain a K x 1 vector in the end               
                channel_over_noise = (abs(bf_weights(:, bb)' * (direct_channel + squeeze(pagemtimes(permute(reflected_channel, [2, 3, 1]), try_conf)))) .^2 / (N0 * B_u)).';       
                obj_prov = V * p_bit * sum(abs(try_conf)) * (1 - alpha) - sum(channel_over_noise .* max(Q_tilde_up, 0));
                % Update configuration
                if obj_prov <= current_obj_prov
                    opt_angle = bb;
                    current_obj_prov = obj_prov;
                end
            end

            %%% looking for the optimal RIS phase-shift            
            % Compute C/N with optimal beamforming
            % transposition to obtain a K x 1 vector in the end            
            channel_over_noise = (abs(bf_weights(:, opt_angle)' * (direct_channel + squeeze(pagemtimes(permute(reflected_channel, [2, 3, 1]), try_conf)))) .^2 / (N0 * B_u)).';
            % Compute greedy obj function
            objective_greedy = V * p_bit * sum(abs(try_conf)) * (1 - alpha) - sum(channel_over_noise .* max(Q_tilde_up, 0));
            % Update configuration
            if objective_greedy <= current_objective_greedy
                current_conf((ll - 1) * N_g + 1 : ll * N_g) = try_conf((ll - 1) * N_g + 1 : ll * N_g);
                current_objective_greedy = objective_greedy;   
                opt_beam = opt_angle;
            end
        end
    end     
end
