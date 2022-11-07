
% Plot the errors from different sims and make convergence plots to check order of accuracy.

close all;
slope1 = 1./NI_list.^1; % 1nd order slope
slope2 = 1./NI_list.^2; % 2nd order slope
slope3 = 1./NI_list.^3; % 3rd order slope

descriptions = ["Cylinder, Dirichlet-Dirichlet", ...
                "Cylinder, Dirichlet-Neumann", ...
                "Cylinder, Neumann-Dirichlet", ...
                "Sphere, Dirichlet-Dirichlet", ...
                "Sphere, Dirichlet-Neumann", ...
                "Sphere, Neumann-Dirichlet"];

for caseNumber = 1:6
    figure(caseNumber);
    slope1plot = loglog(NI_list,slope1*NI_list(1)^1*errors_2Norm(caseNumber,1),'c--', 'LineWidth', 1.);
    hold on
    slope2plot = loglog(NI_list,slope2*NI_list(1)^2*errors_2Norm(caseNumber,1),'b--', 'LineWidth', 1.);
    slope3plot = loglog(NI_list,slope3*NI_list(1)^3*errors_2Norm(caseNumber,1),'g--', 'LineWidth', 1.);
    error1NormPlot = loglog(NI_list,errors_1Norm(caseNumber,:),'ro-', 'LineWidth', 1.5);
    error2NormPlot = loglog(NI_list,errors_2Norm(caseNumber,:),'rx-', 'LineWidth', 1.5);
    xlabel('N');
    ylabel('Error');
    title(descriptions(caseNumber));
    legend([error1NormPlot, error2NormPlot, slope1plot, slope2plot, slope3plot], ["1-norm of error", "2-norm of error", "1st order", "2nd order", "3rd order"]);
end

spreadFigures;