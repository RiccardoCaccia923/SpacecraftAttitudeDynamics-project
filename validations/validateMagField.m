close all
clear 
clc

% set model name
mdl = 'magFieldValidation';  

% define parameters
constant = constantsConfig;
pert = pertConfig;
mu          = constant.mu;          % [m^3/s^2] Terra
Re          = constant.Re;          % [km]
rNorm       = Re + 600;             % [km] radius ~600 km
we          = constant.we;          % [rad/s]
H0          = pert.mf.H0;           % [T]
magAxAngle  = pert.mf.magAxAngle;   % [rad]


% MC settings
N   = 100;                    % # MC cases
rng(1);                       
relErr   = zeros(N,1);
%orthDot  = zeros(N,1);
BRefMC = zeros(N,3);
BSimMC = zeros(N,3);

% open model
load_system(mdl);

for k = 1:N
    % 1) generate random direction (unitary)
    rhat = randn(3,1)';
    rhat = rhat/norm(rhat);

    % 2) define simulation inputs
    in = Simulink.SimulationInput(mdl);
    in = in.setVariable('r_vers', rhat);
    in = in.setVariable('rNorm',  rNorm);
    in = in.setVariable('mu',  mu);
    in = in.setVariable('Re',  Re);
    in = in.setVariable('we',  we);
    in = in.setVariable('H0',  H0);
    in = in.setVariable('magAxAngle',  magAxAngle);

    % 3) istantaneous simulation
    in = in.setModelParameter('StopTime','0');
    simOut = sim(in);

    % 4) read simulated value(To Workspace -> 'BN_sim')
    if isprop(simOut,'BN_sim')
        BN_sim = simOut.BN_sim;
    else
        % fallback: prova dai vari contenitori
        if ~isempty(simOut.find('BN_sim'))
            BN_sim = simOut.find('BN_sim');
        else
            error('cant find ''BN_sim'' in simOut.');
        end
    end
    BN_sim = BN_sim(:).';   % 1x3

    % 5) compute analitic value
    mhat = simOut.mhat;
    BN_ref = H0*(Re/rNorm)^3 * ( 3*dot(mhat,rhat)*rhat - mhat );

    % 6) fill MC parameters
    BRefMC(k,:) = BN_ref;
    BSimMC(k,:) = BN_sim;
    relErr(k) = norm(BN_sim - BN_ref) / max(norm(BN_ref), eps);
    %orthDot(k)= abs(dot(n_B, BN_sim.')) / max(norm(BN_sim), eps);
end

% numeric report
fprintf('Monte Carlo Mag Field (N=%d)\n', N);
fprintf('  RMS relative error   : %.3e\n', sqrt(mean(relErr.^2)));
fprintf('  Max relative error   : %.3e\n', max(relErr));
%fprintf('  < |n · tau_sim| >    : %.3e (should be ~0)\n', mean(orthDot));

% plots
figure; 
semilogy(relErr,'.'); 
grid on;
xlabel('sample #'); 
ylabel('relative error'); 
title('Mag Field Monte Carlo – relative error');

figure; 
scatter3(BRefMC(:,1),BRefMC(:,2),BRefMC(:,3),12,'filled'); 
hold on;
scatter3(BSimMC(:,1),BSimMC(:,2),BSimMC(:,3),12); 
grid on; 
axis equal;
legend('analytic','model'); 
title('B cloud: model vs analytic');
xlabel('B_x'); 
ylabel('B_y'); 
zlabel('B_z');
