function simOptions = simOptions()
    simOptions.SolverType   = 'fixed-step';
    simOptions.Solver       = 'ode5';
    simOptions.FixedStep    = '0.1';
    % simOptions.StartTime    = '0';
    % simOptions.StopTime     = '10';
end