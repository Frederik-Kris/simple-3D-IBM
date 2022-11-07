
% This script defines all input parameters for the code
% Some parameters will be derived from these in setDerivedParameters

caseType = "Cyl"; % Define shell shapes of immersed boundaries. Should be "Cyl" for concentric cylinders or "Sphere" for concentric spheres.
bcType = "DirDir"; % Define BC types at inner and outer shell. Should be "DirDir" for both Dirichlet, "DirNeu" for Neumann on outer, or "NeuDir" for Neumann at inner.
t_end = 1e7; % Stop time, if solution has not converged yet.
convCrit = 1e-10; % Convergence Criterion. 2-norm of change between time levels has to decrease by this factor.
NI  = 41;   % grid resolution in x-direction. The other directions (NJ and NK) are set automatically.
NI_list = [11, 16, 21, 31, 41, 61, 81, 121, 161]; % List of grid resolutions to test if testAllCases is true.
L_x = 1;    % domain size
testAllCases = true; % Overrides caseType and bcType to test all 6 combos. Overrides NI to test increasing resolutions.

innerRadius = 0.2017; % Radii of the shells (cylinders or spheres)
outerRadius = 0.4013;
T_inner = 1; % Temperatures on inner and outer shell. For Dirichlet condition.
T_outer = 2;
tempGradInner = 10; % Temperature gradients on inner and outer shells. For Neumann condition.
tempGradOuter = 10;
xCenter = 0.5012; % Center point for shells. 
yCenter = 0.5011;
zCenter = 0.5013; % z is only used for spheres.
alpha = 1E-6;  % Thermal diffusivity

