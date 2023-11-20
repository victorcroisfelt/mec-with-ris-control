%%% GREEDY RIS OPTIMIZATION FUNCTION %%%
function current_v = Greedy_RIS(current_v, N_blocks, N_RIS, possible_v, K, channel_RIS_up, channel_RIS_AP_up,...
    channel_AP_up, N0, B_u, V, Q_tilde_up, p_bit, alpha)
    N_RIS_block = N_RIS / N_blocks;
    for ll = 1:N_blocks
        current_objective_greedy = +Inf;
        try_v = current_v;
        for vv = 1:length(possible_v)           
            overall_channel_RIS_up = zeros(K, 1); channel_over_noise_up = zeros(K, 1);
            try_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block) = possible_v(vv);
            for kk = 1:K
                overall_channel_RIS_up(kk) = overall_channel_RIS_up(kk) + (channel_RIS_up(:, kk).' * diag(try_v) * channel_RIS_AP_up(:, kk));
                channel_over_noise_up(kk) = abs(channel_AP_up(kk) + overall_channel_RIS_up(kk)) .^2 / (N0 * B_u);
            end

            objective_greedy = V * p_bit * sum(abs(try_v)) * (1 - alpha) - sum(channel_over_noise_up .* max(Q_tilde_up, 0));
            
            if objective_greedy <= current_objective_greedy
                current_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block) = try_v((ll - 1) * N_RIS_block + 1 : ll * N_RIS_block);
                current_objective_greedy = objective_greedy;        
            end
        end
    end     
end