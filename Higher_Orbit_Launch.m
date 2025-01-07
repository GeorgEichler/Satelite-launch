clear all;
close all;

% Define parameters
params.mu = 398600;           % Gravitational parameter for Earth (km^3/s^2)
params.k = .11;               % Turning steepness
params.ro = 6671;             % Orbit radius (300 km above surface)
params.g0 = 0.00981;          % Gravitational acceleration in km/s^2
params.beta = 0.27;           % Contols when turning starts


% Unpack parameters
mu = params.mu;             % Gravitational parameter (G*M)         
k = params.k;               % Turning steepness parameter
ro = params.ro;             % Orbit radius for LEO maneuver
g0 = params.g0;             % Gravitational acceleration (m/s^2)
beta = params.beta;         % Controls when turning starts

% Specific impulse per engine
Isp = [283 421 421 421];                

% Target final orbit
rTarget = 41571;                  

% Delta V Required for Periapsis burn
vPer = sqrt(mu/ro)*(sqrt(2*rTarget/(ro+rTarget))-1);

% Delta V Required for Apoapsis burn
vApo = sqrt(mu/ro)*(1-sqrt((2*ro)/(ro+rTarget)));

% Total delta V for Hohmann transfer
vHoh = vPer + vApo;

% Total delta V for mission
vEnd = sqrt(mu / ro) + vHoh + 3; % target velocity in (km/s)

%structural ratio = structural mass/(structural mass + propellant)
epsilon = [0.0518 0.0808 0.12 0.16]; 

% Mass ratio optimisation
[m0, m, mEmpty, mFuel, mEnd] = Lagrange(vEnd, Isp, epsilon);

% Calculate required thrust for each stage
T = [3*g0*m0 1.8*g0*(m0-m(1)) 1.2*g0*(m0-m(1)-m(2)) 1.2*g0*(m0-m(1)-m(2)-m(3))];

% Calculate mass flow rate of engine for required thrust
mf = (T ./ Isp) / 0.00981;        

%Calculate number of stages 
N = length(epsilon);

% Calculate cumulative mass for each stage
mTotal = zeros(1, length(m));
mTotal(1) = m0; % Total mass at start of the first stage
for i = 2:length(m)
    mTotal(i) = mTotal(i-1) - m(i-1); % Subtract the previous stage mass
end

% Calculate Delta V contributions and check it is enough
DeltaV = (Isp * g0) .* log(mTotal ./ (mTotal - mFuel));
if abs(sum(DeltaV) - vEnd) > 1e-6 % Allow for a small numerical tolerance
    disp('Insufficient Delta V for mission!');
end

% Maneuver 1: Launch to LEO

% Calculate burn required to reach LEO
vLEO = vEnd - vHoh;
vLEO1 = vLEO;
for i = 1:N
    if DeltaV(i) >= vLEO1

        % Solve for mFuel
        mFuelOfFinal = mTotal(i)-mTotal(i)/(exp(vLEO1/(g0*Isp(i))));

        % Store the stage where this occurs
        FinalLEOStage = i;

        % Exit the loop
        break;
    else
        vLEO1 = vLEO1 - DeltaV(i);
    end
end

% Calculate used delta V (Sanity check)
mFuelUpdated = mFuel;
mFuelUpdated(FinalLEOStage) = mFuelOfFinal;
DeltaVUsed = 0;

for i = 1:FinalLEOStage

DeltaVUsed = DeltaVUsed + (Isp(i) * g0) * log(mTotal(i) / (mTotal(i) - mFuelUpdated(i)));

end

% Update masses for spent stages
mFuelUpdated1 = mFuel;
mFuelUpdated1(FinalLEOStage) = mFuel(FinalLEOStage) - mFuelOfFinal;
mTotalUpdated = mTotal;
mTotalUpdated(FinalLEOStage) = mTotal(FinalLEOStage) - mFuelOfFinal;
mUpdated = mTotal;
mUpdated(FinalLEOStage) = m(FinalLEOStage) - mFuelOfFinal;
DeltaVRemaining = 0;

for i = 1 : FinalLEOStage-1
    mTotalUpdated(i) = 0;
    mUpdated(i) = 0;
    mFuelUpdated1(i) = 0;
end

% Calculate remaining delta V (Sanity check)
for i = FinalLEOStage:N

DeltaVRemaining = DeltaVRemaining + (Isp(i) * g0) * log(mTotalUpdated(i) / (mTotalUpdated(i) - mFuelUpdated1(i)));

end


% Calculate remaining delta V for each stage
DeltaVRemaining0 = (Isp * g0) .* log(mTotalUpdated ./ ...
    (mTotalUpdated - mFuelUpdated1));

% Update remaining delta V for empty stages
for i = 1:FinalLEOStage - 1
   DeltaVRemaining0(i)=0; 
end

% Calculate burnout times for LEO launch maneuver
Tb=zeros(1,FinalLEOStage);
for i=1:FinalLEOStage
Tb(1,i)=mFuelUpdated(1,i)/mf(i);
end


% Define initial conditions
r0 = 6371;        % Initial radial position (e.g., 6371 km, Earth's surface)
theta0 = 0;       % Initial angular position (e.g., 0 radians)
vr0 = 0;          % Initial radial velocity (e.g., 0 km/s)
vtheta0 = 0;      % Initial tangential velocity (e.g., 0 km/s)
m_0 = m0;     % Initial mass of the rocket (e.g., 500 kg)

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
for i = 2:FinalLEOStage
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

% Run section and uncomment Plots below to view just the launch to LEO
%Plots
%%

% Now we set orbital contidtions to stay in desired orbit
% This is purely for the sake of accuracy, all variables have to be perfect
% for orbital conditions. A real rocket can self adjust to these perfect
% values, our theoretical rocket cannot so this is sufficient

y0 = y(end, :)';         % Ensure y0 is a column vector
y0(4) = sqrt(mu / ro);   % Set tangential velocity to orbital velocity
y0(3) = 0;               % Set radial velocity to zero
y0(1) = ro;              % Set radial position to Orbit

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

%%
% Maneuver 2: Periapsis Burn

% Entering transfer orbit for first step of Hohmann transfer

% Initial conditions update
y0 = y(end, :)';         % Ensure y0 is a column vector
y0(4) = sqrt(mu / ro) + sqrt(mu/ro)*(sqrt(2*rTarget/(ro+rTarget))-1); % Add
% velocity impulse to attain elliptical transfer orbit
y0(3) = 0;
y0(1) = ro;

% Separate stages spent in Periapsis burn
vPer1 = vPer;
for i = FinalLEOStage:N
    if DeltaVRemaining0(i) <= vPer1
        y0(5) = y0(5) - mUpdated(i);
        vPer1 = vPer1 - DeltaVRemaining0(i);
    else
        % Store when the stage separation occurs
        FinalPerStage = i;
        break;
  
    end
end

% Update masses for spent stages
if FinalPerStage>FinalLEOStage 
for i = FinalLEOStage:FinalPerStage - 1
    mTotalUpdated(i) = 0;
    mUpdated(i) = 0;
    mFuelUpdated1(i) = 0;
end
end

% Calculate remaining delta V for each stage
DeltaVRemaining1 = (Isp * g0) .* log(mTotalUpdated ./ ...
    (mTotalUpdated - mFuelUpdated1));

% Update remaining delta V for empty stages
for i = 1:FinalPerStage-1
   DeltaVRemaining1(i)=0; 
end
DeltaVRemaining1(FinalPerStage) = DeltaVRemaining1(FinalPerStage) + DeltaVRemaining0(FinalPerStage-1) - vPer;


% Time span for simulation and target angle
tspan = [0, 50000];        % Simulate from 0 to t seconds
targetTheta = 3*pi/2;       % Stop at target theta 

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
%%
% Maneuver 3: Apoapsis Burn


% Second velocity impulse to circularise orbit

% Setting initial conditions
y0 = y(end, :)';             % Ensure y0 is a column vector
y0(4) = sqrt(mu / rTarget);  % Setting tangential velocity to orbital 
% velocity of target orbit 
y0(3) = 0;
y0(1) = rTarget;
y0(5) = mEnd;

% Time span for simulation and target angle
tspan = [0, 100000];        % Simulate from 0 to t seconds
targetTheta = 3*pi/2-0.0001;       % Stop at target theta (one complete orbit)

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

   
% Now all_t contains the full time vector and all_y contains the full state trajectory

% Plots the full mission from launch
Plots
