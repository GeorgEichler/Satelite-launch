% Parameters
r_e = 6371;     % Earth's radius in km
r_o = 41571;    % Orbit radius in km
k = 0.0002;     % Growth rate (adjustable) Needs changing depending on orbit, smaller orbit=larger k and vice versa

% Function for f(r) with a delayed transition
f = @(r) (pi/2) * (1 ./ (1 + exp(-k * (r - r_o))));

% Range of r values
r = linspace(r_e, r_o, 1000);

% Trajectory components
x = r .* cos(f(r));
y = r .* sin(f(r));

% Plotting
figure;
plot(x, y, 'b', 'LineWidth', 2);
hold on;

% Adding Earth as a circle
theta = linspace(0, 2*pi, 100);
earth_x = r_e * cos(theta);
earth_y = r_e * sin(theta);
fill(earth_x, earth_y, [0.5, 0.8, 1], 'EdgeColor', 'none'); % Earth fill
plot(earth_x, earth_y, 'k', 'LineWidth', 1); % Earth's edge

% Adding orbit as a circle
orbit_x = r_o * cos(theta);
orbit_y = r_o * sin(theta);
plot(orbit_x, orbit_y, '--r', 'LineWidth', 1); % Orbit

% Labels and visualization settings
axis equal;
xlabel('X (km)');
ylabel('Y (km)');
title('Rocket Trajectory from Launch to Orbit');
legend('Trajectory', 'Earth', 'Orbit');
grid on;

% Debugging: Plot f(r) to visualize the function
figure;
plot(r, f(r), 'g', 'LineWidth', 2);
title('Function f(r)');
xlabel('Radius r (km)');
ylabel('f(r) (radians)');
grid on;

