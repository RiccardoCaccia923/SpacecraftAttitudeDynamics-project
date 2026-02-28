%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR SISO ANALYSIS                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc

Jx = 62;
Jy = 62;
Jz = 81;
%data.Imatrix    = diag([Ix , Iy , Iz]);         % [kg*m^2]

%% Z-AXIS (REACTION WHEEL) -> FAST DYNAMICS
fprintf('\n--- Z-AXIS DESIGN (Reaction Wheel) ---\n');

% A. Requirements
wc_z = 0.04;      % Target Bandwidth [rad/s] (Fast)
PM_z = 80;        % Target Phase Margin [deg] (Robust)
pc_z = 10 * wc_z; % Noise Filter Pole (Far from bandwidth)

% B. Gain Calculation
[Kp_z, Kd_z] = calc_pd(Jz, wc_z, PM_z);
fprintf('Kp = %.6f, Kd = %.6f\n', Kp_z, Kd_z)

% C. Transfer Functions Construction
s = tf('s');
G_z = 1 / (Jz * s^2);                    % Plant (Double Integrator)
C_z = Kp_z + (Kd_z * s) / (1 + s/pc_z);  % Controller (Filtered PD)

% D. Loop Definitions and Labels (TRICK FOR PLOTS)
L_z = C_z * G_z;       % Open Loop
S_z = 1 / (1 + L_z);   % Sensitivity (Error Rejection)
T_z = L_z / (1 + L_z); % Closed Loop (Complementary Sensitivity)
Q_z = C_z / (1 + L_z); 

% Signal Naming
L_z.InputName  = 'e_z';
L_z.OutputName = 'q_z';
S_z.InputName  = 'e_x';
S_z.OutputName = 'q_x';
T_z.InputName  = 'q_{ref}';
T_z.OutputName = 'q_{est}';
Q_z.InputName = 'Noise'; 
Q_z.OutputName = 'Torque Z';

%% X/Y-AXES (MAGNETORQUERS) -> SLOW DYNAMICS

fprintf('\n--- X&Y-AXIS DESIGN (Magnetorquers) ---\n')

% A. Requirements
wc_x = 0.08;      % Target Bandwidth [rad/s]
PM_x = 70;        % Target Phase Margin [deg]
pc_x = 20 * wc_x;

% B. Gain Calculation
[Kp_x, Kd_x] = calc_pd(Jx, wc_x, PM_x);
fprintf('Kp = %.6f, Kd = %.6f\n', Kp_x, Kd_x);

% C. Transfer Functions Construction
G_x = 1 / (Jx * s^2);
C_x = Kp_x + (Kd_x * s) / (1 + s/pc_x);

% D. Loop Definitions and Labels
L_x = C_x * G_x;
S_x = 1 / (1 + L_x);
T_x = L_x / (1 + L_x);
Q_x = C_x / (1 + L_x); 

L_x.InputName  = 'e_x';
L_x.OutputName = 'q_x';
S_x.InputName  = 'e_x';
S_x.OutputName = 'q_x';
T_x.InputName  = 'q_{ref}';
T_x.OutputName = 'q_{est}';
Q_x.InputName = 'Noise'; 
Q_x.OutputName = 'Torque X';

%% PLOTTING AND ANALYSIS
% 1. NYQUIST STABILITY ANALYSIS
figure('Name', 'Nyquist Stability Analysis', 'Color', 'w')
nyquist(L_z, 'b', L_x, 'r')
hold on
grid on
plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2)
legend('L_z (Wheel)', 'L_x (Mag)', 'Critical Point (-1,0)')
%title('Nyquist Plot: Stability Check')
xlabel('Real Axis')
ylabel('Imaginary Axis')
xlim([-2.5 1])
ylim([-1.5 1.5])

% 1. OPEN LOOP COMPARISON (L)
figure('Name', 'Open Loop Comparison (L)', 'Color', 'w')
bode(L_z, 'b', L_x, 'r')
grid on
legend('Z-Axis (Wheel)', 'X-Axis (Mag)')
%title('Open Loop Transfer Function L(s) - Bandwidth Comparison')

% 2. MARGINS
figure('Name', 'Z-Axis Analysis', 'Color', 'w')
margin(L_z,'b')
grid on
hold on
%figure('Name', 'Z-Axis Analysis', 'Color', 'w')
margin(L_x,'r')
%grid on
%title('Stability Margins')
legend('L_z (Wheel)', 'L_x (Mag)')


% 3 SENSITIVITY FUNCTION (S) ANALYSIS
figure('Name', 'Sensitivity S(s) - Robustness')
bodemag(S_z, 'b', S_x, 'r')
grid on
hold on
legend('S_z (Wheel)', 'S_x (Mag)', 'Location', 'SouthEast')
% title('Sensitivity Function S(s)')
Ms_z = getPeakGain(S_z);
Ms_x = getPeakGain(S_x);
% yline(20*log10(Ms_z), 'b--', sprintf('Ms Z: %.2f dB', 20*log10(Ms_z)))
% yline(20*log10(Ms_x), 'r--', sprintf('Ms X: %.2f dB', 20*log10(Ms_x)))
legend('S_z (Wheel)', 'S_x (Mag)', 'Location', 'SouthEast')
fprintf('\n--- SENSITIVITY ANALYSIS (Robustness) ---\n')
fprintf('Peak Sensitivity (Ms) - Target: < 6 dB (2.0)\n')
fprintf('Ms Z-Axis: %.4f (%.2f dB)\n', Ms_z, 20*log10(Ms_z))
fprintf('Ms X-Axis: %.4f (%.2f dB)\n', Ms_x, 20*log10(Ms_x))
if Ms_z > 2.0 || Ms_x > 2.0
    fprintf('WARNING: Peak Sensitivity is high! Reduce Bandwidth or increase Phase Margin.\n')
else
    fprintf('RESULT: Robustness is GOOD (Ms < 2.0).\n')
end

figure('Name', 'Sensitivity S(s) - Performance')
step(S_z, 'b', S_x, 'r')
grid on
hold on
yline(0, 'k--', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left')
legend('Z (Wheel)', 'X (Mag)' , 'Target Reference')
ylabel('Errore Normalizzato')
xlabel('Tempo [s]')

% 4: CLOSED LOOP FUNCTION (COMPLEMENTARY SENSITIVITY FUNCTION T(s))
% Bode
figure('Name', 'Closed Loop Analysis - BODE T(s)')
bodemag(T_z, 'b', T_x, 'r')
grid on
legend('T_z (Wheel)', 'T_x (Mag)', 'Location', 'SouthWest')
% title('Closed Loop Frequency Response (Magnitude)')
yline(-3, 'k--', 'Bandwidth (-3dB)')

% Step Response
figure('Name', 'Closed Loop Analysis - STEP T(s)')
step(T_z, 'b', T_x, 'r')
hold on
grid on
yline(1, 'k--', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left')
% title('Step Response (Closed Loop Performance)')
legend('Response Z (Wheel)', 'Response X (Mag)', 'Target Reference')

fprintf('Closed Loop T(s) Performance')
info_z = stepinfo(T_z);
info_x = stepinfo(T_x);
fprintf('\n Z axis(Wheel):\n')
fprintf(' - Rise Time: %.3f s\n', info_z.RiseTime)
fprintf(' - Overshoot: %.2f %%\n', info_z.Overshoot)
fprintf(' - Bandwidth: %.3f rad/s\n', bandwidth(T_z))

fprintf('\n X axis(Magnetorquer):\n')
fprintf(' - Rise Time: %.3f s\n', info_x.RiseTime)
fprintf(' - Overshoot: %.2f %%\n', info_x.Overshoot)
fprintf(' - Bandwidth: %.3f rad/s\n', bandwidth(T_x))


% 5: CONTROL EFFORT
figure('Name', 'Actuator Effort (Control Sensitivity)', 'Color', 'w');
bodemag(Q_z, 'b', Q_x, 'r');
grid on;
legend('Z (Wheel)', 'X (Mag)');
title('Control Effort (Torque required per radian of error)');

% Step Response

figure('Name', 'Control Effort - Unit Step Response');
step(Q_z, 'b', Q_x, 'r');
grid on;
title('Control Effort 1rad step (1 rad = 57 deg)');
legend('Z (Wheel)', 'X (Mag)');
ylabel('Coppia [Nm]');
xlabel('Tempo [s]');
info_z = stepinfo(Q_z);
info_x = stepinfo(Q_x);
peak_z = info_z.Peak; 
peak_x = info_x.Peak;

fprintf('\n Control Effort (1rad)\n');
fprintf('Peak Z (Wheel): %.4f Nm\n', peak_z);
fprintf('Peak X (Mag):   %.4f Nm\n', peak_x);


%% FUNCTION: GAINs CALCULATION
function [Kp, Kd] = calc_pd(J, wc, PM)
    % Calculates Kp and Kd for a double integrator plant given J, wc, PM
    phase_rad = deg2rad(PM);
    Kd = J * wc;
    Kp = (J * wc^2) / tan(phase_rad);
end