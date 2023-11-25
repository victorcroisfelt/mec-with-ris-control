function [optimal_power_ul, optimal_rate_up, freq_MEH, power_mec, current_conf, power_RIS, opt_beam] = ...
    optimization_RIS(K, B_u, Q_ue_comm, Q_es_proc, Q_virtual, ...
        qtz_RIS, N, N_blocks, channel_ue_ris, channel_ris_ap, channel_ue_ap,...
        N0, V, Pt, possible_f, J, tau_pay, p_bit, kappa, alpha, possible_angles, weights_comb, varargin)

    % Rate Optimization parameters
    Q_tilde_up = B_u * (Q_ue_comm - Q_es_proc + Q_virtual);
    Q_tilde_comp = (Q_es_proc + Q_virtual) .* J;
    m = (0:2^(qtz_RIS)-1).';
    possible_confs = [0; exp(1i*(2*m*pi) / (2^(qtz_RIS)))];

    % Optimize the communication parameters
    current_conf = zeros(N, 1);
    [current_conf, opt_beam] = greedy_optimization(current_conf, N_blocks, N, possible_confs, channel_ue_ris, channel_ris_ap,...
           channel_ue_ap, N0, B_u, V, Q_tilde_up, p_bit, alpha, possible_angles, weights_comb);
        
    if ~isempty(varargin)                        
        current_conf = varargin{1};
    end

    % Compute optimal rate and power
    channel_over_noise_ul = (abs(weights_comb(:, opt_beam)' * (channel_ue_ap + channel_ris_ap * diag(current_conf) * channel_ue_ris)) .^2 / (N0 * B_u)).';
    rate_max_up = min(B_u * log2(1 + channel_over_noise_ul .* Pt), Q_ue_comm / tau_pay);
    p_max = 1 ./ channel_over_noise_ul .* (exp(rate_max_up * log(2) / B_u) - 1);
    optimal_power_ul = min(max(0, Q_tilde_up / (alpha * V * log(2)) - 1 ./ (channel_over_noise_ul)), p_max);          
    % Put to zero power giving negative total normalized queue state
    optimal_power_ul(Q_tilde_up<=0) = 0;
    % Compute optimal rate
    optimal_rate_up = B_u * log2(1 + channel_over_noise_ul .* optimal_power_ul);
    power_RIS = p_bit * sum(abs(current_conf));

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