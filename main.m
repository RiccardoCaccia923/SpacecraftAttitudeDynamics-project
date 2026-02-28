%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT - PROJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc

addpath(genpath(fullfile(pwd,'configs')));
addpath(genpath(fullfile(pwd,'plots')));
addpath(genpath(fullfile(pwd,'miscellaneous')));
addpath(genpath(fullfile(pwd,'validations')));
addpath(genpath(fullfile(pwd,'submodules')));

%% ORBIT DATA
data.a      = 6373+400;
data.e      = 0.2;
data.i      = deg2rad(40);
data.raan   = 0;
data.w      = 0;
data.n      = sqrt(astroConstants(13)/data.a^3);
data.T      = 2*pi/data.n;

%% S/C DATA
data.Ix         = 62;
data.Iy         = 62;
data.Iz         = 81;
data.Imatrix    = diag([data.Ix , data.Iy , data.Iz]);         
data.wL         = [0;0;0];
data.n_sc       = data.n;
data.detumb.wNormMin    = 5e-3;

%% INITIAL CONDITIONS
init.w0      = [8e-2 , 8e-2 , 8e-2];
init.A_BN0   = eye(3);
init.A_LN0   = eye(3);
init.theta   = 0;
init.v0      = 7.5;                      
init.q0      = [0,0,0,1]';
init.lambda0 = 1; 

%% CONFIGURATIONS
simOpt = simOptions();
constant = constantsConfig();
pert = pertConfig();
sensors = sensorsConfig();
MEKF = estimationConfig(sensors,simOpt);
control = controlConfig(data);
actuators = actuatorsConfig();

%% MODEL SIMULATION
out = sim('mainModel.slx'); 

%% POST PROCESS
% conservation(data.Imatrix,out.tout,out.w)
plots(out,data)

