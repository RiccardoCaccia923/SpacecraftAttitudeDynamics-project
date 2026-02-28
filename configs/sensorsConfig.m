function sensors = sensorsConfig()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                CONFIGURATION FILE FOR SENSORS PARAMETERS                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% magnetometer
    sensors.magnetometer.frequency          = 10;                                   % [Hz]
    sensors.magnetometer.natFreq            = 2*pi*sensors.magnetometer.frequency;  % [Hz]
    sensors.magnetometer.sampleTime         = 1/sensors.magnetometer.frequency;     % [s]
    sensors.magnetometer.delay              = 1;                                    % [#step]
    sensors.magnetometer.bias               = [50,50,50]*1e-9;                      % 50[nT]
    sensors.magnetometer.SFlinear           = 1667e-6;                              % 5000[ppm] (3sigma)
    sensors.magnetometer.SFnonLinear        = 166.7e-6;                             % 500[ppm] (3sigma)
    sensors.magnetometer.SFnonSimm          = 16.7e-6;                              % 50[ppm] (3sigma)
    sensors.magnetometer.noiseDensity       = 10e-9;                                % 10[nT]
    sensors.magnetometer.Rmisalign          = misalignMat(deg2rad(0.33));           % 1[°] (3sigma)
    sensors.magnetometer.MnonOrthog         = nonOrtMat(deg2rad(0.33));             % 1[°] (3sigma)
    sensors.magnetometer.saturation         = 60000e-9;                             % ±60000[nT] 
    sensors.magnetometer.resolution         = 7.324e-9;                             % 7.324[nT]
    sensors.magnetometer.alpha              = 3e-5;

    sensors.magnetometer.noiseSeed          = [1 2 3];

% star sensor
    sensors.starSensor.frequency            = 5;                                    % [Hz]
    sensors.starSensor.sampleTime           = 1/sensors.starSensor.frequency;       % [s]
    sensors.starSensor.angNoise             = 4.85e-5;                              % 30[arcsec] (3sigma) [rad]=''×π/(180×3600)
    sensors.starSensor.noiseDensity         = 1e-4;
    sensors.starSensor.MsensFrame           = misalignMat(deg2rad(0.1));
    sensors.starSensor.boresightAxis        = [0,0,1]';                             % z axis boresight
    sensors.starSensor.FOV                  = deg2rad(60);
    sensors.starSensor.maxMag               = 6;
    sensors.starSensor.alpha                = 0.25;

    sensors.starSensor.noiseSeed            = [1 2 3];
    
    [sensors.starSensor.starVecs_N ,sensors.starSensor.MAG] = ...
    starsFromCatalogue('starCat50.tsv',sensors.starSensor.maxMag);
end
