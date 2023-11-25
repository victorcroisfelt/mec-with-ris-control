function [R] = eval_rate(channel_ue_ris, channel_ris_ap, channel_ue_ap, N0, B_u, configuration, ap_combiner, ue_power)
%%% Compute the channel capacity given the channels, the power, the 
% configuration and the combiner loaded

channel_over_noise = (abs(ap_combiner' * (channel_ue_ap + channel_ris_ap * diag(configuration) * channel_ue_ris)) .^2 / (N0 * B_u)).';
R = B_u * log2(1 + channel_over_noise .* ue_power);