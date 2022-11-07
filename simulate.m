
% Runs ONE simulation, with ONE set o parameters. "main" calls this multiple times for multiple simulations.

setDerivedParameters;

% Initialize solution array:
T = 1.5*ones(NI,NJ,NK)+0;    % T = Temperature

setNodeCategoryFlags; % Flagging and finding the indices to different node categories, (fluid, solid, boundary, etc)

findGhostAndImagePoints; % Find the ghost nodes among the solid nodes, and find their respective image points

t  = 0;           % t = time
intermSolutn = T; % Intermediate solution, which we will compute using slopes k1,k2,k3.
normHistory = 1;  % Will be an array with the 2-norm of change between each time level. To check convergence.
timeLevel = 0;    % Time level. Increases by one each time step.
tic;              % Start timer
reportInterval=4; % How often to write status report to screen. Period in seconds.
nextReportTime = reportInterval;
previousProgress = 0; % Used to estimate remaining runtime in 'statusReport'

while t+dt <= t_end && normHistory(end) > normHistory(1) * convCrit % Until we hit end-time or convergance criterion
    marchTimeStep;
    if toc > nextReportTime % If it's time to report to screen, do it.
        statusReport;
    end
end
if t < t_end && normHistory(end) > normHistory(1) * convCrit % If not converged, do an adapted time step to hit end-time exactly.
    dt = t_end - t;
    marchTimeStep;
end
T(:,:,:) = NaN; % I want NaN on the solid nodes, except ghost nodes.
updateSurfaceNodeValues; % Update surface- and ghost node values, to match the final time level
updateGhostNodeValues;
T([fluidNodesIndices, ghostNodesIndices]) = intermSolutn([fluidNodesIndices, ghostNodesIndices]); % Write ghost and fluid node values back to the main array.
statusReport; % One last status report, to give total wall-clock time, and say if it converged.
computeAnalyticalSolution;