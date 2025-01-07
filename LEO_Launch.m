clear all;
close all;

% Define parameters
params.mu = 398600;           % Gravitational parameter for Earth (km^3/s^2)
params.k = .36 ;              % Sigmoid steepness
params.ro = 6671;             % Orbit radius (300 km above surface)
params.g0 = 0.00981;          % Gravitational acceleration in km/s^2
params.beta = 0.15;            % Controls when turning starts


% Unpack parameters
mu = params.mu;             % Gravitational parameter (G*M)         
k = params.k;               % Turning steepness parameter
ro = params.ro;             % Orbit radius for LEO maneuver
g0 = params.g0;             % Gravitational acceleration (m/s^2)
beta = params.beta;         % Controls when turning starts

% Specific impulse per engine
Isp = [283 421];                  

% Target final orbit
%rTarget = 41571;                  

% target velocity in (km/s)
vEnd = sqrt(mu / ro) + 3; 


%structural ratio = structural mass/(structural mass + propellant)
epsilon = [0.0518 0.0808]; 


% Lagrange mass ratio optimisation
[m0, m, mEmpty, mFuel, mEnd] = Lagrange(vEnd, Isp, epsilon);

% Calculate required thrust for each stage
T = [3*g0*m0 1.2*g0*(m0-m(1))];

% Mass flow rate of engine
mf = (T ./ Isp) / 0.00981;        

%Calculate number of stages 
N = length(epsilon);

% Calculate Delta V contribution of each stage
DeltaV = (Isp * g0) .* log((m) ./ (mEmpty));
sum(DeltaV);

%Calculate burnout times
Tb=zeros(1,N);
for i=1:N
Tb(1,i)=mFuel(1,i)/mf(i);
end

% Define initial conditions
r0 = 6371;        % Initial radial position (e.g., 6371 km, Earth's surface)
theta0 = 0;       % Initial angular position (e.g., 0 radians)
vr0 = 0;          % Initial radial velocity (e.g., 0 km/s)
vtheta0 = 0;      % Initial tangential velocity (e.g., 0 km/s)
m_0 = m0;         % Initial mass of the rocket (e.g., 500 kg)

% Pack initial conditions into a column vector
y0 = [r0; theta0; vr0; vtheta0; m_0];

% Time span for simulation of first stage
tspan = [0, Tb(1,1)];        % Simulate from 0 to t seconds

% Initialize storage for results
all_t = []; % To store time points
all_y = []; % To store state variables


% First stage simulation
mfN = mf(1);
IspN = Isp(1);
[t, y] = ode15s(@(t, y) Equations_Motion(t, y, params,IspN,mfN), [0, Tb(1,1)], y0);

% Append results
all_t = [all_t; t];
all_y = [all_y; y];


% Loop for subsequent stages
for i = 2:N
    % Update initial conditions for the next stage
    y0 = y(end, :)'; % Ensure y0 is a column vector
    y0(5) = y0(5) - mEmpty(i-1);

    
    % Time span for simulation of this stage
    tspan = [0, Tb(1, i)];
    
    % Solve the ODE system
    mfN = mf(i);
    IspN = Isp(i);
    [t, y] = ode15s(@(t, y) Equations_Motion(t, y, params,IspN,mfN), tspan, y0);
    
    % Append results
    all_t = [all_t; all_t(end) + t]; % Shift time to continue from previous stage
    all_y = [all_y; y];
end

%% 

% Now we set orbital contidtions to stay in desired orbit
% This is purely for the sake of accuracy, all variables have to be perfect
% for orbital conditions. A real rocket can self adjust to these perfect
% values, our theoretical rocket cannot so this is sufficient

% Initial Conditions 
y0 = y(end, :)';         % Ensure y0 is a column vector
y0(4) = sqrt(mu / ro);
y0(3) = 0;
y0(1) = ro;
y0(5) = mEnd;

% Time span for simulation and target angle
tspan = [0, 4000];        % Simulate from 0 to t seconds
targetTheta = pi/2;       % Stop at target theta 

% Pass parameters and events to ode45
options = odeset('Events', @(t, y) thetaEvent(t, y, params, targetTheta));

% Solve the ODE system
[t, y, te, ye, ie] = ode45(@(t, y) Orbital_Equations(t, y, params), tspan, y0, options);

% Output: 
% t    -> Time points
% y    -> Solution values
% te   -> Time when the event occurred
% ye   -> State values when the event occurred
% ie   -> Index of the triggered event

% Append results
    all_t = [all_t; all_t(end) + t]; % Shift time to continue from previous stage
    all_y = [all_y; y];

% Generate Plots
Plots