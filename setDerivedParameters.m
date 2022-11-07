
% Derives some parameters based on the input ones. Plus some sanity check.

if (caseType ~= "Cyl" && caseType ~= "Sphere")
    disp("Invalid case type. Must be ""Cyl"" or ""Sphere"".");
end
if (bcType ~= "DirDir" && bcType ~= "DirNeu" && bcType~= "NeuDir")
    disp("Invalid BC type. Must be ""DirDir"", ""DirNeu"" or ""NeuDir"".");
end
NJ  = NI;
L_y = L_x;
if caseType == "Cyl"
    NK  = 3;
    L_z = 2/(NI-1);
elseif caseType == "Sphere"
    NK = NI;
    L_z = L_x;
end
epsilon = eps(L_x) + eps(L_y) + eps(L_z); % Buffer for machine precision.
dx = L_x / (NI-1);  % Grid spacings
dy = L_y / (NJ-1);
dz = L_z / (NK-1);
if (dx ~= dy || dx ~= dz)
    disp("Choose grid and domain size to give dx=dy=dz");
end
dt = dx^2 / (6*alpha); % Stability limit if dx=dy=dz