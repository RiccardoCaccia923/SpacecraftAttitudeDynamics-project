function actuators = actuatorsConfig()

    actuators.magTorq.maxDipole = 30;                % [Am^2]

    actuators.wheel.maxTorque   = 0.265;             % [Nm]
    actuators.wheel.maxMomentum = 18;                % [Nms]
    actuators.wheel.maxSpeed    = 4000 * (2*pi/60);  % 4000 [RPM]

end

