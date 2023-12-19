clear all; clc;

%% Common values
M = 8;
N = 64;
K = 4;

save_dir = 'data/plots/';

setup = randi(10);
M_sample = 50;

% Take one to load all the variables
file_name = ['data/results/ce/N', num2str(N), '_M', num2str(M), '_K', num2str(K), '_N_block32_discount0.1.mat'];
load(file_name);
num_block = length(N_blocks_vec);
num_discount = length(discount_vec);

%% Loading data
L_one = zeros(num_block, num_discount, num_arrival, N_slot);
E_one = zeros(num_block, num_discount, num_arrival, N_slot);

L_avg = zeros(num_block, num_discount, num_arrival, num_points);
E_avg = zeros(num_block, num_discount, num_arrival, num_points);

ind_conf_one = 1;

ind_v = lyapunov_tradeoff == 1e8;

% Blocks
for ind_block = 1:num_block
    % Discount factor
    for ind_dis = 1:num_discount
        % Load file
        file_name = ['data/results/ce/N', num2str(N), '_M', num2str(M), '_K', num2str(K), '_N_block', num2str(N_blocks_vec(ind_block)), '_discount', num2str(discount_vec(ind_dis), '%1.1f'), '.mat'];
        load(file_name);        
        
        % Arrival rate
        for ind_arr = 1:num_arrival    
            % Load one shot results
            L_one(ind_block, ind_dis, ind_arr, :) = squeeze(max(avg_delay(ind_arr, setup, ind_conf_one, :, :, ind_v), [],  4));
            E_one(ind_block, ind_dis, ind_arr, :) = squeeze(mean(total_energy(ind_arr, setup, ind_conf_one, :, ind_v), 2));
            
            % Load average results
            L_avg(ind_block, ind_dis, ind_arr, :) = squeeze(mean(max(avg_delay(ind_arr, :,:,:,end-M_sample:end, ind_v), [], 4), [2,5])) ;
            E_avg(ind_block, ind_dis, ind_arr, :) = squeeze(mean(total_energy(ind_arr, :,:,end-M_sample:end, ind_v), [2, 4]));
            lab_arr{ind_arr} = ['$\bar{A}/\tau = $ ', num2str(avg_arrival_rate(ind_arr)/1000), ' kbit/s'];
        end 

        lab_dis{ind_dis} = ['$\mu = $ ', num2str(discount_vec(ind_dis), '%1.1f')];
    end
end

%% One shot (one per block and data rate)
close all;

for ind_block = 1
    for ind_arr = 3
        % Generate figures
        figure;
        sgtitle(['$\bar{A}/\tau = $ ', num2str(avg_arrival_rate(ind_arr)/1000), ' kbit/s, $N_g = $', num2str(N/N_blocks_vec(ind_block))], 'interpreter', 'latex');
       
        % energy
        subplot(2, 1, 1);        
        plot(squeeze(E_one(ind_block, :, ind_arr, :)).');                
        grid on;
        xlabel('slot $t$', 'interpreter', 'latex')
        ylabel('Total energy $E_\sigma^\mathrm{tot}(t)$ [J]', 'interpreter', 'latex');
        legend(lab_dis, 'interpreter', 'latex');

        % Latency
        subplot(2, 1, 2);   
        hold on;
        plot(squeeze(L_one(ind_block, :, ind_arr, :)).');                
        hold off;    
        grid on;
        xlabel('slot $t$', 'interpreter', 'latex')
        ylabel('$\max_{k\in\mathcal{K}} L_k(t)$ [ms]', 'interpreter', 'latex');
        legend(lab_dis, 'interpreter', 'latex');
    end
end
filename = [save_dir, 'ce_oneshot_', 'N', num2str(N), '_M', num2str(M), '_K', num2str(K)];
savefig(filename);
matlab2tikz([filename, '.tex']);

%% Total
close all;

ind_dis = 2; 

figure(99);
sgtitle('Overall performance CE', 'interpreter', 'latex')
subplot(2, 2, 1);
title('$N_g = 2$', 'interpreter', 'latex')
plot(conf_codebook_size, squeeze(E_avg(1, :, 2, :))); 
grid on;
xlabel('$C_\mathrm{ce}$', 'interpreter', 'latex')        
ylabel('$E_\sigma^\mathrm{tot}(t)$ [J]', 'interpreter', 'latex'); 
legend(lab_dis, 'Interpreter','latex');
  
subplot(2, 2, 3);
plot(conf_codebook_size, squeeze(L_avg(1, :, 2, :)) / 1000);
grid on;
xlabel('$C_\mathrm{ce}$', 'interpreter', 'latex')        
ylabel('$\max_{k\in\mathcal{K}} L_k(t)$ [s]', 'interpreter', 'latex'); 
set(gca,'yscale','log')
line([conf_codebook_size(1), conf_codebook_size(end)], 0.3 *[1, 1], 'color', 'black', 'LineStyle', '--')


subplot(2, 2, 2);
title('$N_g = 1$', 'interpreter', 'latex')
plot(conf_codebook_size, squeeze(E_avg(2, ind_dis, :, :))); 
grid on;
xlabel('$C_\mathrm{ce}$', 'interpreter', 'latex')        
ylabel('$E_\sigma^\mathrm{tot}(t)$ [J]', 'interpreter', 'latex'); 
  
subplot(2, 2, 4);
plot(conf_codebook_size, squeeze(L_avg(2, ind_dis, :, :)) / 1000);
grid on;
xlabel('$C_\mathrm{ce}$', 'interpreter', 'latex')        
ylabel('$\max_{k\in\mathcal{K}} L_k(t)$ [s]', 'interpreter', 'latex'); 
set(gca,'yscale','log')
line([conf_codebook_size(1), conf_codebook_size(end)], 0.3*[1, 1], 'color', 'black', 'LineStyle', '--')

filename = [save_dir, 'ce_', 'N', num2str(N), '_M', num2str(M), '_K', num2str(K)];
savefig(filename);
matlab2tikz([filename, '.tex']);