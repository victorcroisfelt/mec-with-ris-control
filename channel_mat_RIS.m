function [Hdir, H1, H2] = channel_mat_RIS( ...
        Nt, Nr, N_RIS, ...
        dist_ap_ris, angl_ap_ris, ...
        dist_minimum, dist_maximum, ...
        f, no_mat, varargin)

    %% Parameters
    lambda = 3e8/f;     % Wavelength
    dt = lambda/2;      % TX antenna space
    dr = lambda/2;      % RX antenna space
    dris = lambda/2;    % RIS element space
    k = 2*pi/lambda;    % Wavenumber
    kappa = 0;          % Rician kappa constant !!! CHECK LATER

    % Load the FSPL of the direct link
    if isempty(varargin)                        
        alpha = 2;
    else
        alpha = varargin{1};
    end 

    %% Geometrical placement, x, y and z axis 

    % AP antenna array
    rx_arr(1,:) = ones(1, Nr);
    rx_arr(2,:) = (sort(0:Nr-1,'descend')-(Nr-1)/2) * dr; 
    rx_arr(3,:) = zeros(1, Nr);

    % Get AP position
    x_pos_ap = dist_ap_ris * -sin(deg2rad(angl_ap_ris));
    y_pos_ap = dist_ap_ris * cos(deg2rad(angl_ap_ris));

    rx_arr(1,:) = x_pos_ap * rx_arr(1,:);
    rx_arr(2,:) = rx_arr(2,:) + y_pos_ap;    

    % UE antenna array
    tx_arr(1,:) = zeros(1, Nt); 
    tx_arr(2,:) = (sort(0:Nt-1, 'descend')-(Nt-1)/2) * dt; 
    tx_arr(3,:) = zeros(1,Nt); 

    % Generate UE position
    distance = sqrt(rand() * (dist_maximum^2 - dist_minimum^2) + dist_minimum^2);
    angle = pi/2 * rand();

    x_pos_ue = distance * sin(angle);
    y_pos_ue = distance * cos(angle);
    
    tx_arr(1,:) = tx_arr(1,:) + x_pos_ue; 
    tx_arr(2,:) = tx_arr(2,:) + y_pos_ue; 

    % RIS 
    center = [0 0]; % RIS center position 
    N1 = sqrt(N_RIS);
    N2 = N1;                                    % Number of RIS elements in two dimensions N1 and N2
    ris_pos = RISPosition(N1,N2,dris,center);   % RIS elements' coordinates 
    a = repmat(ris_pos{1},N1,1);                

    % Placing RIS elements in proper coordinates
    ris_arr(1,:) = a(:)';        
    ris_arr(2,:) = zeros(1,N_RIS);
    ris_arr(3,:) = repmat(ris_pos{2},1,N2); 

    %% Compute Distances

    % Distance between the TX and RX antennas
    for i1 = 1:Nr                                                                                   
        for j1 = 1:Nt
            dist_ap_ue_antennas(i1,j1) = norm(rx_arr(:,i1)-tx_arr(:,j1));
        end 
    end 

    % Indirect paths (TX-RIS-RX)
    for l1 = 1:N_RIS  % Distance between the RIS elements and the RX antennas                                                            
        for r1 = 1:Nr  
            dist_ap_ris_antennas(r1,l1) = norm(rx_arr(:,r1)-ris_arr(:,l1)); 
        end
        for t1 = 1:Nt % Distance between the RIS elements and the TX antennas                  
            dist_ue_ris_antennas(l1,t1) = norm(tx_arr(:,t1)-ris_arr(:,l1));   
        end
    end

    %% Get Channels
    
    % Direct link, LOS matrix exponents
    Hdir_los = exp(-1i * k * dist_ap_ue_antennas);   
    
    % Inversion of the FSPL of the direct link 
    FSPL_dir = (lambda/(4*pi))^2 / (sqrt(distance^2 + dist_ap_ris^2) ^ alpha(1));    

    % Direct link channel matrix
    Hdir = Rician(Hdir_los, sqrt(FSPL_dir), no_mat, kappa);              
    
    % TX-RIS channel matrix
    H1_los = exp(-1i * k * dist_ue_ris_antennas);
    FSPL_ind1 = (lambda/(4*pi))^2 / (distance ^ 2);
    FSPL_1 = sqrt(FSPL_ind1);
    H1 = Rician(H1_los, FSPL_1, no_mat, kappa);

    % RIS-RX channel matrix
    H2_los = exp(-1i * k * dist_ap_ris_antennas); % RIS-RX link, LOS matrix exponents 
    FSPL_ind2 = (lambda/(4*pi))^2 / (dist_ap_ris ^ 2);
    FSPL_2 = sqrt(FSPL_ind2);
    H2 = Rician(H2_los, FSPL_2, no_mat, kappa);
end