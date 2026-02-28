function [Kp, Kd] = PDcontrollerGains(J, wc, PM)
    % Calculates Kp and Kd for a double integrator plant given J, wc, PM
    phase_rad = deg2rad(PM);
    Kd = J * wc;
    Kp = (J * wc^2) / tan(phase_rad);
end