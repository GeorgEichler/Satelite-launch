function dydt = Orbital_Equations(t, y, params)
    % Unpack state variables
    r = y(1);         % Radial position
    theta = y(2);     % Angular position
    vr = y(3);        % Radial velocity
    vtheta = y(4);    % Tangential velocity
    m = y(5);         % Mass of the rocket
    
    % Unpack parameters
    mu = params.mu;             % Gravitational parameter (G*M)         
    k = params.k;               % Turning steepness parameter
    ro = params.ro;             % Orbit radius
    g0 = params.g0;             % Gravitational acceleration (m/s^2)
    beta = params.beta;         % Controls when turning starts
   

 
    % Define ODEs
    drdt = vr;                              % Radial position
    dthetadt = vtheta / r;                  % Angular position
    dvrdt = -(mu / r^2) + vtheta^2/r;       % Radial velocity
    dvthetadt = - vr * vtheta / r;       % Tangential acceleration

    % Mass update
    dmdt = 0;  % Mass
    
    % Return derivatives as a column vector
    dydt = [drdt; dthetadt; dvrdt; dvthetadt; dmdt];
end

