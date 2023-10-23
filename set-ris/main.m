clear all; clc; rand('state',0); randn('state',0)

% Define distance from UE to the RIS
D = 500;

% Vector of probabilities
proba_vec = linspace(0, 1, 11);

% Vector of angles
angles_vec = (pi/2) * linspace(0, 1, 11);

% Prepare to save simulation results
avg_delay = zeros(length(proba_vec), length(angles_vec), 4);
rate = zeros(length(proba_vec), length(angles_vec), 4);

% Simulation
for pp = 1:length(proba_vec)
    proba = proba_vec(pp);
    
    tic
    for aa = 1:length(angles_vec)
        angle = angles_vec(aa);

        [avg_delay(pp, aa, :), rate(pp, aa, :)] = RIS_MEC_Control_UL_siso(D, angle, proba);

    end
    elapsed_time = toc;

    disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);

end

string = ['data/set-ris_D' + num2str(D) + '.mat'];

save(string)