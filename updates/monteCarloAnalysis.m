%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MONTE CARLO CONFIGURATION
MC.n_runs = 50;                 
MC.enable_plots = false;         
MC.results_APE = zeros(MC.n_runs, 1);
MC.results_RPE = zeros(MC.n_runs, 1);
MC.settling_time = zeros(MC.n_runs, 1);

fprintf('Starting Monte Carlo Analysis(%d runs)...\n', MC.n_runs)

%% MONTE CARLO LOOP
for i = 1:MC.n_runs
    
    fprintf('Run %d/%d... \n', i, MC.n_runs)
    
    % 1.a Initial condition - Angular Velocity
    % Esempio: Variazione random +/- 20% sulla w0 nominale
    perturbation_w = (rand(1,3) - 0.5) * 0.4; % fattore 0.4 crea range [-0.2, +0.2]
    init.w0 = nominal_w0 .* (1 + perturbation_w);
    
    % 1.b Initial Condition - Attitude
    rand_axis = rand(3,1) - 0.5;
    rand_axis = rand_axis / norm(rand_axis);
    max_error_deg = 10; 
    rand_angle_rad = deg2rad(max_error_deg * (2*rand - 1));
    K = [0 -rand_axis(3) rand_axis(2); rand_axis(3) 0 -rand_axis(1); -rand_axis(2) rand_axis(1) 0];
    C_noise = eye(3) + sin(rand_angle_rad)*K + (1-cos(rand_angle_rad))*(K*K);
    init.A_BN0 = C_noise * eye(3);

    % 1.c Initial Condition - quaternion
    % Genera un quaternione randomico uniforme
    % rand_axis = rand(3,1) - 0.5;
    % rand_axis = rand_axis / norm(rand_axis);
    % rand_angle = rand * 2 * pi; % Rotazione casuale 0-360 deg
    % q_rand = [sin(rand_angle/2)*rand_axis; cos(rand_angle/2)]; % esempio internet 
    % init.q0 = q_rand; 
    
    % 2. Inertia Uncertainty
    J_unc = (rand(3,1) - 0.5) * 0.2; 
    data.Ix = nominal_Imatrix(1,1) * (1 + J_unc(1));
    data.Iy = nominal_Imatrix(2,2) * (1 + J_unc(2));
    data.Iz = nominal_Imatrix(3,3) * (1 + J_unc(3));
    data.Imatrix = diag([data.Ix, data.Iy, data.Iz]);
    
    % 3. Sensor Noise
    sensors.magnetometer.noiseSeed = randi([1, 100], 1, 3);
    sensors.starSensor.noiseSeed = randi([1, 100], 1, 3);
    
    % Simulation
    out = sim('mainModel.slx', 'SrcWorkspace', 'current');
    
    % Results
    t_end = out.tout(end);
    idx_steady = out.tout > (t_end - 2000); % Ultimi 2000 secondi
    
    % APE
    MC.results_APE(i) = mean(out.thetaErr(idx_steady));

    % RPE
    MC.results_RPE(i) = rad2deg(mean(out.w_norm(idx_steady)));

    % Attitude Estimation Error
    MC.results_AEE(i) = mean(out.thetaErrMekf(idx_steady));

    % Settling Time Detumbling 
    threshold = data.detumb.wNormMin;
    idx_settled = find(out.w_norm > threshold, 1, 'last');
    if isempty(idx_settled)
        MC.settling_time(i) = 0;
    else
        MC.settling_time(i) = out.tout(idx_settled);
    end

end

disp('Analisi Monte Carlo Completata.')

%% MC PLOTS & STATISTICS
figure('Name', 'Monte Carlo Analysis', 'Color', 'w')

% 1. Detumbling Settling Time
subplot(2,2,1)
histogram(MC.settling_time/data.T, 15, 'FaceColor', [0.2 0.4 0.8])
xline(5, 'r--', 'LineWidth', 1.5)
title('Detumbling Settling Time')
xlabel('T orb')
ylabel('Occurrences')
grid on

% 2. APE (Accuracy) 
subplot(2,2,2)
histogram(MC.results_APE, 20, 'FaceColor', [0.8 0.4 0.2])
xline(5, 'r--', 'LineWidth', 1.5) % Linea Requisito APE
title('Steady-State APE (Accuracy)')
xlabel('Error [deg]')
grid on

% 3. RPE (Stability) 
subplot(2,2,3)
histogram(MC.results_RPE, 20, 'FaceColor', [0.2 0.7 0.3])
xline(0.1, 'r--', 'LineWidth', 1.5) % Linea Requisito RPE
title('Steady-State RPE (Stability)')
xlabel('Error [deg/s]')
grid on

% 4. Attitude Estimation Error (AEE)
subplot(2,2,4)
histogram(MC.results_AEE, 20, 'FaceColor', [0.2 0.7 0.3])
xline(1, 'r--', 'LineWidth', 1.5) % Linea Requisito AEE
title('Attitude Determination Error')
xlabel('Error [deg]')
grid on

% 4. (APE vs RPE)
figure
scatter(MC.results_APE, MC.results_RPE, 20, 'filled')
yline(0.1, 'r--')
xline(5, 'r--') % Box dei requisiti
title('Performance Correlation')
xlabel('APE [deg]')
ylabel('RPE [deg/s]')
grid on

fprintf('\n OVERALL RESULTS\n');
fprintf('Mean Settling Time: %.2f s\n', mean(MC.settling_time));
fprintf('Mean APE:           %.4f deg (Max: %.4f)\n', mean(MC.results_APE), max(MC.results_APE));
fprintf('Mean RPE:           %.4f deg/s (Max: %.4f)\n', mean(MC.results_RPE), max(MC.results_RPE));
fprintf('Mean AEE:           %.4f deg/s (Max: %.4f)\n', mean(MC.results_AEE), max(MC.results_AEE));


% Verifica Pass/Fail
pass_APE = sum(MC.results_APE < 5);
pass_RPE = sum(MC.results_RPE < 0.1);
pass_AEE = sum(MC.results_AEE < 1);

fprintf('Compliance APE:     %d/%d runs (%.1f%%)\n', pass_APE, MC.n_runs, (pass_APE/MC.n_runs)*100);
fprintf('Compliance RPE:     %d/%d runs (%.1f%%)\n', pass_RPE, MC.n_runs, (pass_RPE/MC.n_runs)*100);
fprintf('Compliance AEE:     %d/%d runs (%.1f%%)\n', pass_AEE, MC.n_runs, (pass_AEE/MC.n_runs)*100);