function dydt = Equations_Motion(t, y, params,IspN,mfN)
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

    %Calculate Thrust (KgKm/s^2)
    T=IspN*mfN*g0;

    % Midpoint for transition
    r_mid = 6371 + beta * (ro - 6371);

    % Orbital velocity
    v_orbit = sqrt(mu / ro);

    % Turning Function
    alpha = (pi/2) / (1 + exp(-k * (r - r_mid)));  % Smooth transition to orbit
    
    % Define ODEs
    drdt = vr;  % Radial position
    dthetadt = vtheta / r;  % Angular position
    dvrdt = -mu / r^2 + (T / m) * cos(alpha) + vtheta^2/r;  % Radial velocity
    dvthetadt = (T / m) * sin(alpha)  - vr * vtheta / r;  % Tangential velocity

    % Mass update
    dmdt = -mfN;  % Mass
    
    % Return derivatives as a column vector
    dydt = [drdt; dthetadt; dvrdt; dvthetadt; dmdt];
end

