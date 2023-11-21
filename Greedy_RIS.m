%%% GREEDY RIS OPTIMIZATION FUNCTION %%%
function [current_v, opt_beam] = Greedy_RIS( ...
    current_v, N_blocks, N_RIS, possible_v, K, channel_RIS, channel_RIS_AP,...
    channel_AP, N0, B_u, V, Q_tilde_up, p_bit, alpha, possible_angles, ...
    weights_comb, Nr)

    N_RIS_block = N_RIS / N_blocks;
    for ll = 1:N_blocks
        current_objective_greedy = +Inf; current_obj_prov = +Inf;
        try_v = current_v;
        for vv = 1:length(possible_v)           
            overall_channel_RIS = zeros(Nr, K); channel_over_noise = zeros(K, 1);
            try_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block) = possible_v(vv);

            for bb = 1:possible_angles
                for kk = 1:K
                    overall_channel_RIS(:, kk) = overall_channel_RIS(:, kk) + (channel_RIS_AP(:, :, kk) * diag(try_v) * channel_RIS(:, kk));
                    channel_over_noise(kk) = abs(weights_comb(:, bb)' * (channel_AP(:, kk) + overall_channel_RIS(:, kk))) .^2 / (N0 * B_u);
                end
                obj_prov = V * p_bit * sum(abs(try_v)) * (1 - alpha) - sum(channel_over_noise .* max(Q_tilde_up, 0));
                % looking for the optimal AP beamforming
                if obj_prov <= current_obj_prov
                    opt_angle = bb;
                    current_obj_prov = obj_prov;
                end
            end

            for kkk = 1:K
                overall_channel_RIS(:, kkk) = overall_channel_RIS(:, kkk) + (channel_RIS_AP(:, :, kkk) * diag(try_v) * channel_RIS(:, kkk));
                channel_over_noise(kkk) = abs(weights_comb(:, opt_angle)' * (channel_AP(:, kkk) + overall_channel_RIS(:, kkk))) .^2 / (N0 * B_u);
            end
            objective_greedy = V * p_bit * sum(abs(try_v)) * (1 - alpha) - sum(channel_over_noise .* max(Q_tilde_up, 0));
            % looking for the optimal RIS phase-shift
            if objective_greedy <= current_objective_greedy
                current_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block) = try_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block);
                current_objective_greedy = objective_greedy;   
                opt_beam = opt_angle;
            end
        end
    end     
end
