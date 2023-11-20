clear all; clc; rand('state',0); randn('state',0)

% Define distance from UE to the RIS
D = 500;

% Vector of angles
angles_vec = (pi/2) * linspace(0, 1, 11);

% Prepare to save simulation results
avg_delay = zeros(length(angles_vec), 4);
rate = zeros(length(angles_vec), 4);

% Simulation
%for pp = 1:length(proba_vec)
%    proba = proba_vec(pp);
    

tau = 400e-3
k = 4;
T = k * 1/14 * 10e-3;

tic
for aa = 1:length(angles_vec)
    angle = angles_vec(aa);

    [avg_delay(aa, :), rate(aa, :)] = RIS_MEC_Control_UL_siso(D, angle, tau, T);

end
elapsed_time = toc;
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);

%end

string = ['data/baseline.mat'];
save(string)