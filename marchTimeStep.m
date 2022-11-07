% Everything that needs to happen to bring the solution from time t to t+dt

% Step 1:
updateSurfaceNodeValues;
updateGhostNodeValues;
k1 = RK4_step(intermSolutn, alpha, activeNodesIndices, NI, NJ, dx, dy, dz); % 1st RK4 step. Finding the slope at time t using Euler's method
intermSolutn(activeNodesIndices) = T(activeNodesIndices) + dt/2 * k1; % Find solution at t+dt/2 using slope k1.
domainBC; % Apply BC at the boundaries of the domain, NOT the immersed boundary.
% Step 2:
updateSurfaceNodeValues;
updateGhostNodeValues;
k2 = RK4_step(intermSolutn, alpha, activeNodesIndices, NI, NJ, dx, dy, dz); % 2nd RK4 step. Finding the slope at time t+dt/2
intermSolutn(activeNodesIndices) = T(activeNodesIndices) + dt/2 * k2; % Find solution at t+dt/2 using slope k2.
domainBC; % Apply BC at the boundaries of the domain, NOT the immersed boundary.
% Step 3:
updateSurfaceNodeValues;
updateGhostNodeValues;
k3 = RK4_step(intermSolutn, alpha, activeNodesIndices, NI, NJ, dx, dy, dz); % 3rd RK4 step. Finding the slope at time t+dt/2 again
intermSolutn(activeNodesIndices) = T(activeNodesIndices) + dt * k3; % Find solution at t+dt using slope k3.
domainBC; % Apply BC at the boundaries of the domain, NOT the immersed boundary.
% Step 4:
updateSurfaceNodeValues;
updateGhostNodeValues;
k4 = RK4_step(intermSolutn, alpha, activeNodesIndices, NI, NJ, dx, dy, dz); % 4th RK4 step. Finding the slope at time t+dt
intermSolutn(activeNodesIndices) = T(activeNodesIndices) + dt/6 * (k1 + 2*k2 + 2*k3 + k4); % Find solution at t+dt using all the slopes k1,...,k4.
domainBC; % Apply BC at the boundaries of the domain, NOT the immersed boundary.
normHistory(timeLevel+1) = sqrt(sum((intermSolutn(fluidNodesIndices)-T(fluidNodesIndices)).^2)); % Compute the 2-norm of the difference between T(t,...) and T(t+dt,...), and add it to the history. It is 
T(fluidNodesIndices) = intermSolutn(fluidNodesIndices); % Write solution back to the main temperature array

timeLevel = timeLevel + 1;
t = t + dt;


