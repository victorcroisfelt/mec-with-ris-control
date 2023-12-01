clear all; close all; clc;

%% Common values
M = 8;
N = 64;
K = 2;

save_dir = 'data/plots/';

setup = randi(10);
M_sample = 70;


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
    L = squeeze(mean(avg_delay(1, setup, 1:4, :, :, 1), 4));
    E = squeeze(total_energy(1, setup, 1:4, :, 1)) * 1000;
    
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
    L = squeeze(mean(max(avg_delay(1, :,:,:,end-M_sample:end, 1), [], 4), [2, 5]));
    E = squeeze(mean(total_energy(1, :,:,end-M_sample:end, 1), [2, 4])) * 1000;

    figure(5);
    subplot(2, 1, 1);
    if error_code == 0
        hold on;   
        line([prob_vector(2), prob_vector(end)], [E(1) E(1)], 'color', 'black', 'linestyle', '--');
        grid on;
        xlabel('error probability $p_k^i$', 'interpreter', 'latex')
        ylabel('Total Energy $E_\sigma^\mathrm{tot}$ [mJ]', 'interpreter', 'latex');        
    end    
    plot(prob_vector(2:end), smoothdata(E(2:end))); 
    if error_code == 3        
        legend('error-free', 'INI-U', 'INI-R', 'SET-U', 'SET-R', 'interpreter', 'latex');
        set(gca,'xscale','log')
    end

    subplot(2, 1, 2);
    if error_code == 0        
        hold on;     
        line([prob_vector(2), prob_vector(end)], [L(1) L(1)], 'color', 'black', 'linestyle', '--');        
        grid on;
        xlabel('error probability $p_k^i$', 'interpreter', 'latex')        
        ylabel('Avg. latency [ms]', 'interpreter', 'latex');
    end
    plot(prob_vector(2:end), L(2:end)); 
    
    if error_code == 3        
        legend('error-free', 'INI-U', 'INI-R', 'SET-U', 'SET-R', 'interpreter', 'latex');
        set(gca,'xscale','log')
        set(gca,'yscale','log')        
    end    
end

figure(5); hold off; 
sgtitle('Overall performance error proba', 'interpreter', 'latex')


filename = [save_dir, 'error_proba_', 'N', num2str(N), '_M', num2str(M), '_K', num2str(K)];
savefig(filename);
matlab2tikz([filename, '.tex']);