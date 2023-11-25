%%% GREEDY OPTIMIZATION FUNCTION for RIS and COMBINER %%%
function [current_conf, opt_beam] = greedy_optimization( ...
    current_conf, N_blocks, N, possible_confs, ...
    channel_ue_ris, channel_ris_ap, channel_ue_ap, ...
    N0, B_u, V, Q_tilde_up, p_bit, alpha, possible_angles, ...
    weights_comb)

    N_RIS_block = N / N_blocks;

    % Checking for each block of RIS elements
    for ll = 1:N_blocks
        current_objective_greedy = +Inf; current_obj_prov = +Inf;
        try_conf = current_conf;
        for vv = 1:length(possible_confs)           
            try_conf((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block) = possible_confs(vv);
            
            %%% looking for the optimal AP combiner
            for bb = 1:possible_angles
                % Compute C/N with combiner and configuration under test
                % transposition to obtain a K x 1 vector in the end
                channel_over_noise = (abs(weights_comb(:, bb)' * (channel_ue_ap + channel_ris_ap * diag(try_conf) * channel_ue_ris)) .^2 / (N0 * B_u)).';                
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
            channel_over_noise = (abs(weights_comb(:, opt_angle)' * (channel_ue_ap + channel_ris_ap * diag(try_conf) * channel_ue_ris)) .^2 / (N0 * B_u)).';
            % Compute greedy obj function
            objective_greedy = V * p_bit * sum(abs(try_conf)) * (1 - alpha) - sum(channel_over_noise .* max(Q_tilde_up, 0));
            % Update configuration
            if objective_greedy <= current_objective_greedy
                current_conf((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block) = try_conf((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block);
                current_objective_greedy = objective_greedy;   
                opt_beam = opt_angle;
            end
        end
    end     
end
