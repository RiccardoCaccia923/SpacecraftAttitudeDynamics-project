close all
clear
clc

% set model name
mdl = 'ggPertValidation'; 

% defne parameters
constant    = constantsConfig;
pert        = pertConfig;
mu          = constant.mu;          % [m^3/s^2] Terra
Re          = constant.Re;          % [km]
rNorm       = Re + 600;             % [km] radius ~600 km
J           = diag([0.4 0.6 0.8]);           % [kg m^2] inertia in body
K           = 3*mu/rNorm^3;                      % teorical scale factor

% MC settings
N   = 100;                    % # MC cases
rng(1);                       
relErr   = zeros(N,1);
orthDot  = zeros(N,1);
tauRefMC = zeros(N,3);
tauSimMC = zeros(N,3);

% open model
load_system(mdl);

for k = 1:N
    % 1) generate casual nadir direction (unitary) in BODY
    n_B = randn(3,1); n_B = n_B/norm(n_B);

    % 2) define simulation inputs
    in = Simulink.SimulationInput(mdl);
    in = in.setVariable('n_B', n_B);
    in = in.setVariable('J',   J);
    in = in.setVariable('mu',  mu);
    in = in.setVariable('Re',  Re);
    in = in.setVariable('rNorm',   rNorm);

    % 3) istantaneous simulation
    in = in.setModelParameter('StopTime','0');
    simOut = sim(in);

    % 4) read simulated values
    if isprop(simOut,'tau_sim')
        tau_sim = simOut.tau_sim;
    else
        if ~isempty(simOut.find('tau_sim'))
            tau_sim = simOut.find('tau_sim');
        else
            error('cant find ''tau_sim'' in simOut');
        end
    end
    tau_sim = tau_sim(:).';   % 1x3

    % 5) compute analytic value
    tau_ref = (K * cross(n_B, J*n_B)).';

    % 6) fill MC parameters
    tauRefMC(k,:) = tau_ref;
    tauSimMC(k,:) = tau_sim;
    relErr(k) = norm(tau_sim - tau_ref) / max(norm(tau_ref), eps);
    orthDot(k)= abs(dot(n_B, tau_sim.')) / max(norm(tau_sim), eps);
end

% numeric report
fprintf('Monte Carlo GG (N=%d)\n', N);
fprintf('  RMS relative error   : %.3e\n', sqrt(mean(relErr.^2)));
fprintf('  Max relative error   : %.3e\n', max(relErr));
fprintf('  < |n · tau_sim| >    : %.3e (should be ~0)\n', mean(orthDot));

% Plots
figure
semilogy(relErr,'.') 
grid on
xlabel('sample #')
ylabel('relative error')
title('GG Monte Carlo – relative error')

figure
scatter3(tauRefMC(:,1),tauRefMC(:,2),tauRefMC(:,3),12,'filled') 
hold on
scatter3(tauSimMC(:,1),tauSimMC(:,2),tauSimMC(:,3),12)
grid on
axis equal
legend('analytic','model')
title('GG torque cloud: model vs analytic')
xlabel('\tau_x') 
ylabel('\tau_y') 
zlabel('\tau_z')
