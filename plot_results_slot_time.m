clear all; close all; clc;

%% Common values
M = 8;
N = 64;
K = 4;

save_dir = 'data/plots/';

setup = randi(10);
M_sample = 50;

%% Plotting time-slot
% Load file
file_name = ['data/results/slot_time/N', num2str(N), '_M', num2str(M), '_K', num2str(K), '_last.mat'];
load(file_name);
num_setups_mean = 1:num_setups;

lab = cell(num_lyapunov, 1);
for n=1:num_lyapunov
    lab{n} = ['$V$ = ', num2str(lyapunov_tradeoff(n), '%1.0e')];
end

for ind_arr = 1:num_arrival
    for ind_tau = [1, 3, 11]
        % Take the metrics
        L_one = squeeze(max(avg_delay(ind_arr, setup, ind_tau, :, :, :), [],  4));
        E_one = squeeze(mean(total_energy(ind_arr, setup, ind_tau, :, :), 2)) * 1000;
        
        % Generate figure
        figure;
        sgtitle(['$AR_\mathrm{avg} = $ ', num2str(avg_arrival_rate(ind_arr)/1000), ' kbit/s, $\tau = $ ', num2str(tau_range(ind_tau))], 'interpreter', 'latex');
        subplot(2, 1, 1);
        hold on;
        for n=1:num_lyapunov
            plot(E_one(:, n));
        end
        hold off;
        grid on;
        xlabel('slot $t$', 'interpreter', 'latex')
        ylabel('Total energy $E_\sigma^\mathrm{tot}(t)$ [mJ]', 'interpreter', 'latex');
        legend(lab, 'interpreter', 'latex');
        
        subplot(2, 1, 2);   
        hold on;
        for n=1:num_lyapunov
            plot(L_one(:, n)); 
        end
        hold off;    
        grid on;
        xlabel('slot $t$', 'interpreter', 'latex')
        ylabel('Max. latency [ms]', 'interpreter', 'latex');
        legend(lab, 'interpreter', 'latex');

    end

ind_v = lyapunov_tradeoff == 1e8;

% Load average results
L = squeeze(mean(max(avg_delay(ind_arr, num_setups_mean,:,:,end-M_sample:end, ind_v), [], 4), [2,5])) ;
E = squeeze(mean(total_energy(ind_arr, num_setups_mean,:,end-M_sample:end, ind_v), [2, 4])) * 1000 ;

figure(99);
subplot(3, 1, 1);
plot(tau_range * 1000, smoothdata(E) / 1000); hold on;

subplot(3, 1, 2);
plot(tau_range * 1000, smoothdata(E) ./ tau_range.' / 1000); hold on;
   
subplot(3, 1, 3);
plot(tau_range * 1000, smoothdata(L)); hold on;

label{ind_arr} = ['$AR_\mathrm{avg} = $', num2str(avg_arrival_rate(ind_arr) / 1000), ' kbits/s'];

end 

figure(99);
sgtitle('Overall performance slot time', 'interpreter', 'latex')
subplot(3, 1, 1);       
grid on;
xlabel('slot duration $\tau$ [ms]', 'interpreter', 'latex')        
ylabel('$E_\sigma^\mathrm{tot}(t)$ [J]', 'interpreter', 'latex'); 
legend(label, 'Interpreter','latex');

subplot(3, 1, 2);       
grid on;
xlabel('slot duration $\tau$ [ms]', 'interpreter', 'latex')        
ylabel('$E_\sigma^\mathrm{tot}(t)/\tau$ [W]', 'interpreter', 'latex'); 
legend(label, 'Interpreter','latex');


subplot(3, 1, 3);
line(tau_range * 1000, 300 * ones(length(tau_range)), 'color', 'black', 'LineStyle', '--')
grid on;
xlabel('slot duration $\tau$ [ms]', 'interpreter', 'latex')        
ylabel('$\max_k L(t)$ [ms]', 'interpreter', 'latex'); 

figure(99);
filename = [save_dir, 'time_slot_', 'N', num2str(N), '_M', num2str(M), '_K', num2str(K)];
savefig(filename);
matlab2tikz([filename, '.tex']);