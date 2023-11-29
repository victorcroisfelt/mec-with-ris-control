clear all; close all; clc;

%% Common values
M = 8;
N = 64;
K = 4;

save_dir = 'data/plots/';

setup = randi(100);
colors = {'blue', 'red', 'magenta', 'green'};
M_sample = 50;

%% Plotting error proba results
for error_code = 0:3

    % Take filename
    switch error_code
        case 0
            error_type = 'ini-u';                    
        case 1
            error_type = 'ini-r';                    
        case 2
            error_type = 'set-u';                    
        case 3
            error_type = 'set-r';                    
     end
    
    % Load file
    file_name = ['data/results/', error_type, '/N', num2str(N), '_M', num2str(M), '_K', num2str(K), '.mat'];
    load(file_name);

    %%% check results for a single setup with perfect an unperfect control    

    L = squeeze(mean(avg_delay(4, setup, 1:4, :, :), 4));
    E = squeeze(total_energy(4, setup, 1:4, :)) * 1000;
    
    figure(error_code + 1);
    sgtitle(['Metrics vs iteration for ' error_type, ' setup n. ', num2str(setup)], 'interpreter', 'latex')

    subplot(2, 1, 1);
    hold on;
    for n=1:4
        plot(E(n, :));
    end
    hold off;
    grid on;
    xlabel('slot $t$', 'interpreter', 'latex')
    ylabel('Total energy $E_\sigma^\mathrm{tot}$ [mJ]', 'interpreter', 'latex');
    legend('$p_k^{i}$ = 0', '$p_k^{i}$ = 0.001', '$p_k^{i}$ = 0.01', '$p_k^{i}$ = 0.1', 'interpreter', 'latex');

    subplot(2, 1, 2);   
    hold on;
    for n=1:4
        plot(L(n, :)); 
    end
    hold off;    
    grid on;
    xlabel('slot $t$', 'interpreter', 'latex')
    ylabel('Avg. latency [ms]', 'interpreter', 'latex');
    legend('$p_k^{i}$ = 0', '$p_k^{i}$ = 0.001', '$p_k^{i}$ = 0.01', '$p_k^{i}$ = 0.1', 'interpreter', 'latex');

    % savefig([save_dir, 'oneshot_', 'N', num2str(N), '_M', num2str(M), '_K', num2str(K), '_', error_type]);

    %%% Load average results
    L = squeeze(mean(avg_delay(4, :,:,:,end-M_sample:end), [2, 4, 5]));
    E = squeeze(mean(total_energy(4, :,:,end-M_sample:end), [2, 4])) * 1000;

    figure(5);
    subplot(2, 1, 1);
    if error_code == 0
        hold on;   
        line(prob_vector(2:end), repmat(E(1), length(prob_vector(2:end)), 1), 'color', 'black', 'linestyle', '--');
        grid on;
        xlabel('error probability $p_k^i$', 'interpreter', 'latex')
        ylabel('Total Energy $E_\sigma^\mathrm{tot}$ [mJ]', 'interpreter', 'latex');        
    end    
    plot(prob_vector(2:end), E(2:end)); 
    if error_code == 3        
        legend('error-free', 'INI-U', 'INI-R', 'SET-U', 'SET-R', 'interpreter', 'latex');
        set(gca,'xscale','log')
    end

    subplot(2, 1, 2);
    if error_code == 0        
        hold on;     
        line(prob_vector(2:end), repmat(L(1), length(prob_vector(2:end)), 1), 'color', 'black', 'linestyle', '--');
        grid on;
        xlabel('error probability $p_k^i$', 'interpreter', 'latex')        
        ylabel('Avg. latency [ms]', 'interpreter', 'latex');
    end
    plot(prob_vector(2:end), L(2:end)); 
    
    if error_code == 3        
        legend('error-free', 'INI-U', 'INI-R', 'SET-U', 'SET-R', 'interpreter', 'latex');
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        ylim([200 3000]);
    end    
end

figure(5); hold off; 
sgtitle('Overall performance', 'interpreter', 'latex')
savefig([save_dir, 'error_proba_', 'N', num2str(N), '_M', num2str(M), '_K', num2str(K)]);


%% Plotting time-slot
close all;

% Load file
file_name = ['data/results/slot_time/N', num2str(N), '_M', num2str(M), '_K', num2str(K), '_error_type0_try2.mat'];
load(file_name);

ind_tau = tau_range == 0.05 | tau_range == 0.15 | tau_range == 0.3;
% ind_arr = A_avg_vector == 256;

for ind_arr = 1:5

L_one = squeeze(mean(avg_delay(ind_arr, setup, ind_tau, :, :), 4));
E_one = squeeze(mean(total_energy(ind_arr, setup, ind_tau, :), 2)) * 1000;

figure;
sgtitle(['$A_\mathrm{avg} = $ ', num2str(A_avg_vector(ind_arr))], 'interpreter', 'latex');
subplot(2, 1, 1);
hold on;
for n=1:3
    plot(E_one(n, :));
end
hold off;
grid on;
xlabel('slot $t$', 'interpreter', 'latex')
ylabel('Total energy $E_\sigma^\mathrm{tot}$ [mJ]', 'interpreter', 'latex');
legend('$\tau$ = 50 [ms]', '$\tau$ = 150 [ms]', '$\tau$ = 300 [ms]', 'interpreter', 'latex');

subplot(2, 1, 2);   
hold on;
for n=1:3
    plot(L_one(n, :)); 
end
hold off;    
grid on;
xlabel('slot $t$', 'interpreter', 'latex')
ylabel('Avg. latency [ms]', 'interpreter', 'latex');
legend('$\tau$ = 50 [ms]', '$\tau$ = 150 [ms]', '$\tau$ = 300 [ms]', 'interpreter', 'latex');



% Load average results
L = squeeze(mean(max(avg_delay(ind_arr, :,:,:,end-M_sample:end), [], 4), [2,5]));
E = squeeze(mean(total_energy(ind_arr, :,:,end-M_sample:end), [2, 4])) * 1000;

figure(10);
subplot(2, 1, 1);
plot(tau_range * 1000, E); hold on;
   
subplot(2, 1, 2);
plot(tau_range * 1000, L); hold on;

label{ind_arr} = ['$A_\mathrm{avg} = $', num2str(A_avg_vector(ind_arr))];

end 

figure(10);
sgtitle('Overall performance', 'interpreter', 'latex')
subplot(2, 1, 1);       
grid on;
xlabel('slot duration $\tau$ [ms]', 'interpreter', 'latex')        
ylabel('Total Energy $E_\sigma^\mathrm{tot}$ [mJ]', 'interpreter', 'latex'); 
legend(label, 'Interpreter','latex');

subplot(2, 1, 2);
plot(tau_range * 1000, tau_range * 2500, 'color', 'black', 'LineStyle', '--')
grid on;
xlabel('slot duration $\tau$ [ms]', 'interpreter', 'latex')        
ylabel('Avg. latency ', 'interpreter', 'latex'); 

savefig([save_dir, 'time_slot_', 'N', num2str(N), '_M', num2str(M), '_K', num2str(K)]);