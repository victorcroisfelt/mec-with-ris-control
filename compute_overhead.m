function [tau_overhead, perc, tau_sig, tau_ce, tau_ra] = compute_overhead(tti_time, tau_ris, tau, N_pilots, N_blocks, conf_codebook_size, qtz_RIS, M, f_ra, K, ris_cc)
% Function computing the overhead time depending on:
% tti_time, time of pilots and control packets (T in the paper);
% tau_ris, the switching time of the RIS (\tau_s in the paper);
% tau, total duration of the slot:
% N_pilots, number of pilots per user (N_p in the paper);
% N_blocks, number of RIS block to optimize (N/N_g in the paper);
% qtz_RIS, number of quantization bits per RIS element (b in the paper);
% M, number of AP's antennae
% f_ra, frequency allocated for RA
% K, number of UEs;


% Compute signaling time
if strcmp(ris_cc, 'ib-cc')
    tau_sig = 3 * tau_ris + 5 * tti_time; 
else
    tau_sig = 3 * tau_ris + 3 * tti_time; 
end

% Compute channel estimation overhead
tau_ce = (tau_ris + tti_time * N_pilots) * (conf_codebook_size + 1);

% Compute mu constant
mu = 2 * (K * (3*N_blocks + M + 4) + 3);

% Compute RA time 
n_ra = (2.^qtz_RIS + 1) * N_blocks * mu;
tau_ra =  n_ra / f_ra;

% Compute overhead time
tau_overhead = tau_sig + tau_ce + tau_ra;

if any(tau <= tau_overhead)
    error('at least one tau value lower than tau_overhead')
end

% Percentage of slot duration dedicated to offloading
perc = 1 - (tau_overhead ./ tau);

end

