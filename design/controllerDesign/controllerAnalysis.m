%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR SISO ANALYSIS                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc

%% 1. define system parameters
% I = 0.08; %data.Imatrix(3,3);
I = [0.04,0.06,0.08];
data.Imatrix = diag(I);
wc = 10;
zc = 0.1;
pc = 100;
%%
for  j = 1:3

    I = data.Imatrix(j,j);

    % |L| = k * |j*wc + zc| / ( I * |j*wc|^2 * |j*wc + pc| ) = 1
    num = I * abs(1i*wc)^2 * abs(1i*wc + pc);
    den = abs(1i*wc + zc);
    k = num/den;

    % 2. define transfer functions
    s = tf('s');
    
    C(j) = k * (s + zc) / (s + pc);    % controller C(s) = k * (s + zc) / (s + pc)
    
    G(j) = 1 / (I * s^2);              % plant G(s) = 1 / (I*s^2)
    
    L(j) = C(j) * G(j);                      % Loop T.F. L(s) = C(S) * G(s)
    
    S(j) = 1 / (1 + L(j));                % Sensitivity Function
    
    F(j) = L(j) / (1 + L(j));                % Complementary Sensitivity Function
    
    Q(j) = C(j) / (1 + L(j));                % Control Sensitivity Function
end

%% 3. stability analysis
% --- FIGURA 1: NYQUIST (Stabilità) ---
figure('Name', 'Nyquist Plot');
nyquist(L(1)); 
hold on
nyquist(L(2));
nyquist(L(3));
grid on;
title('Nyquist Plot di L(s)');
xlim([-2 2]);
ylim([-2 2]);
%%
% --- FIGURA 2: BODE (Margini) ---
figure('Name', 'Bode & Margini X')
margin(L(1))
grid on

figure('Name', 'Bode & Margini Y')
margin(L(2))
grid on

figure('Name', 'Bode & Margini Z')
margin(L(3))
grid on

%%

figure
bodemag(F)
grid on

% --- FIGURA 3: CONTROL EFFORT (Sforzo Attuatori) ---
figure('Name', 'Control Effort (KS)');
bodemag(Q, 'r'); 
grid on;
title('Control Effort (Coppia richiesta per Rumore unitario)');
% ylabel('Magnitudine [Nm / rad]');
% xlabel('Frequenza [rad/s]');

%% step response
figure
step(S(1))

figure 
step(F(1))