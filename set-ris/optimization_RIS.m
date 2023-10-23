%%% OPTIMIZATION %%%
function [optimal_power_up, optimal_rate_up, freq_MEH, power_mec, current_v, power_RIS] = ...
    optimization_RIS( ...
        K, B_u, Q_local, Q_MEH, Z, qtz_RIS, N_RIS, N_blocks, channel_RIS_up,...
        channel_RIS_AP_up, channel_AP_up, N0, V, Pt, possible_f, J, delta,...
        p_bit, kappa, alpha, proba)

    % Rate Optimization
    Q_tilde_up = B_u * (Q_local - Q_MEH + Z);
    Q_tilde_comp = (Q_MEH + Z) .* J;
    m = [0:2^(qtz_RIS)-1]';
    possible_v = [0; exp(1i*(2*m*pi) / (2^(qtz_RIS)))];

    current_v = zeros(N_RIS, 1);
    current_v = Greedy_RIS(current_v, N_blocks, N_RIS, possible_v, K, channel_RIS_up, channel_RIS_AP_up,...
           channel_AP_up, N0, B_u, V, Q_tilde_up, p_bit, alpha);
    
    proba = proba;
    toss_coins = rand(K, 1);
    config_v = zeros(N_RIS, 1);%dftmtx();

    overall_channel_RIS_up = zeros(K, 1); channel_over_noise_up = zeros(K, 1);
    for kk = 1:K
        
        if toss_coins(kk) < proba      
            overall_channel_RIS_up(kk) = overall_channel_RIS_up(kk) + (channel_RIS_up(:, kk).' * diag(config_v) * channel_RIS_AP_up(:, kk));
        else
            overall_channel_RIS_up(kk) = overall_channel_RIS_up(kk) + (channel_RIS_up(:, kk).' * diag(current_v) * channel_RIS_AP_up(:, kk));
        end
        channel_over_noise_up(kk, 1) = abs(channel_AP_up(kk) + overall_channel_RIS_up(kk)) .^2 / (N0*B_u);
        
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