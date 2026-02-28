%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR SISO ANALYSIS                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc

Ix = 0.04;
Iy = 0.06;
Iz = 0.08;
data.Imatrix    = diag([Ix , Iy , Iz]);         % [kg*m^2]

%% Define System Parameters
J_mat = data.Imatrix; 
Jx = J_mat(1,1);
Jy = J_mat(2,2);
Jz = J_mat(3,3);

% Banda Passante (wc): Quanto veloce deve reagire?
wc_wheel = 0.1;    % [rad/s] Asse Z (Veloce, attuatore forte)
wc_mag   = 0.02;   % [rad/s] Assi X/Y (Lenti, attuatori deboli)

% Margine di Fase (PM): Quanto smorzato? (60-70 deg è l'ottimo)
PM_req   = 70;     % [deg]

%% Compute Gain Z axis (Wheel) - (Yaw)
fprintf('\n--- ASSE Z (Reaction Wheel) ---\n');

% Compute Phase Factor (equal for all axis)
phase_factor = tan(deg2rad(PM_req)); 

% Derivative: (Fissa la banda passante: Kd*wc / J*wc^2 = 1)
Kd_z = Jz * wc_wheel;

% Proportional: (fix phase margin moving the zero)
Kp_z = (Jz * wc_wheel^2) / phase_factor;

fprintf('Jz: %.4f kgm^2 | Banda: %.2f rad/s\n', Jz, wc_wheel);
fprintf('Kp_z = %.6f\n', Kp_z);
fprintf('Kd_z = %.6f\n', Kd_z);

%% Compute gain X/Y axis (MagTorquers) - (Roll/Pitch)
fprintf('\n--- ASSI X/Y (Magnetorquers) ---\n');

% X axis (Roll)
Kd_x = Jx * wc_mag;
Kp_x = (Jx * wc_mag^2) / phase_factor;

fprintf('Jx: %.4f | Kp_x = %.6f | Kd_x = %.6f\n', Jx, Kp_x, Kd_x);

% Y axis (Pitch)
Kd_y = Jy * wc_mag;
Kp_y = (Jy * wc_mag^2) / phase_factor;

fprintf('Jy: %.4f | Kp_y = %.6f | Kd_y = %.6f\n', Jy, Kp_y, Kd_y);

%% Check for Saturation
err_test_rad = deg2rad(10);

torque_req_z = Kp_z * err_test_rad;
torque_req_x = Kp_x * err_test_rad;

fprintf('\n--- CHECK SATURATION (static error 10 deg) ---\n');
fprintf('Requsted Torque on Z (Wheel): %.4f Nm\n', torque_req_z);
fprintf('Requsted Torque on X/Y (Mag):   %.4f Nm\n', torque_req_x);

if torque_req_x > 1e-3 % Magnetorquers tipical limit
    warning('Attention: Magnetic Gain is Too High, reduce wc_mag.');
else
    disp('Magneic check: OK (Probably it doesnt saturate instantly)');
end
