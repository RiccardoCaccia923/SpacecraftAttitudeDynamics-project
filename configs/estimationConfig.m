function MEKF = estimationConfig(sensors,simOptions)

    MEKF.ssNoiseVariance = sensors.starSensor.noiseDensity/sqrt(2*sensors.starSensor.sampleTime)*1e-2;
    MEKF.magNoiseVariance = sensors.magnetometer.noiseDensity/sqrt(2*sensors.magnetometer.sampleTime)*1e5;
    MEKF.qTheta = 1e-10;
    MEKF.qOmega = 1e-7;
    MEKF.noiseParams = [MEKF.ssNoiseVariance ; MEKF.magNoiseVariance ; MEKF.qTheta ; MEKF.qOmega];
    MEKF.init.P = eye(6)*1e-4;
    MEKF.init.q = [1,0,0,0];
    MEKF.init.w = [0,0,0];
    MEKF.dt = str2double(simOptions.FixedStep);

    MEKF.F = [zeros(3) , eye(3) ; 
              zeros(3) , zeros(3)];

    MEKF.PHI = eye(6) + MEKF.F * MEKF.dt;

    MEKF.Q = [eye(3,3)*MEKF.qTheta , zeros(3,3) ; 
              zeros(3,3) , eye(3,3)*MEKF.qOmega] * MEKF.dt;

    R_ss_mat = [1,0,0;0,1,0;0,0,4] * MEKF.ssNoiseVariance; 
    R_mag_mat = eye(3) * MEKF.magNoiseVariance;
    MEKF.R = blkdiag(R_ss_mat, R_mag_mat);

end