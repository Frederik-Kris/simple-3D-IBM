% --------------------------------------------------------------
% |  3D Immersed Boundary Method test
% --------------------------------------------------------------

clear;

setParameters; % Settings to define the test case(s).

if testAllCases
    errors_1Norm = zeros(6, length(NI_list)); % Arrays to store 1 and 2-norms of errors for different sims. Magic constant 6 is the number of combos of caseType and bcType (2x3).
    errors_2Norm = zeros(6, length(NI_list));
    currentCaseNumber = 1;
    for caseType = ["Cyl", "Sphere"] % Test all combos of caseType and bcType
        for bcType = ["DirDir", "DirNeu", "NeuDir"]
            currentGridNumber = 1;
            for NI = NI_list % Test all resolutions in NI_list
                simulate;
                currentGridNumber = currentGridNumber + 1;
            end
            currentCaseNumber = currentCaseNumber + 1;
        end
    end
    makeConvergencePlots; % Don't plot all solutions from the simulations, only convergence plots.
else
    simulate;
    makeFiguresSingleSim; % Plot solution and error (and comparison with Martin's code if cylinder).
end







