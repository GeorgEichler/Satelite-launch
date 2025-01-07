function [value, isterminal, direction] = thetaEvent(t, y, params, targetTheta)
% Event function for stopping at target theta

    theta = y(2);  % Extract the angular position
    value = theta - targetTheta; % Event when theta reaches targetTheta
    isterminal = 1;             % Stop the integration
    direction = 0;              % Detect both increasing and decreasing
end