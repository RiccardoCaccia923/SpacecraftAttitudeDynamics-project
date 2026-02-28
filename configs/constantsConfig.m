function constants = constantsConfig()

    constants.G       = astroConstants(1);             % [km^3/(kg*s^2)]
    constants.Me      = 5.972*1e24;                    % [kg]
    constants.Re      = astroConstants(23);            % [km]
    constants.mu      = astroConstants(13);            % [km^3/s^2]
    constants.we      = 7.29*1e-5;                     % [rad/s]
    constants.rSun    = astroConstants(2);
    constants.radSun  = astroConstants(81);
    constants.c       = astroConstants(5)*1e3;         % [m/s]

end