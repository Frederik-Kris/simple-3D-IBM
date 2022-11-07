% problem parameters
sigma = -1.0;    % heat production source term: Laplace(T) = sigma
Tinner = 12;     % temperature on inner circle for Dirichlet boundary
Touter = 6;      % temperature on outer circle for Dirichlet boundary
qNeumann = -2.53; % dT/dn on inner or outer Neumann boundary. positive: heat flow into domain, negative: heat flow out of domain
doPlot = false;

% storage & formatting
nMin = 40;
nMax = 250;
nPoints = 32;
nValues = floor(logspace(log10(nMin),log10(nMax),nPoints));
eValues = zeros(size(nValues));
eValues_1storder = zeros(size(nValues));
eValues_2ndorder = zeros(size(nValues));
leg{1} = 'IBM';
leg{2} = 'First Order';
leg{3} = 'Second Order';
sty{1} = 'r-';
sty{2} = 'k--';
sty{3} = 'b--';

figure(1);

results = zeros(nPoints, 5);

for iCase = 1:3

    if iCase == 1
        BCTypeInner = 'Dirichlet';
        BCTypeOuter = 'Dirichlet';
    elseif iCase == 2
        BCTypeInner = 'Neumann';
        BCTypeOuter = 'Dirichlet';
    else
        BCTypeInner = 'Dirichlet';
        BCTypeOuter = 'Neumann';
    end
    

    for i = 1:length(nValues)
        eValues(i) = laplace(nValues(i), BCTypeOuter, BCTypeInner, sigma, Tinner, Touter, qNeumann, doPlot);
        results(i, iCase) = eValues(i);
        if (iCase == 1)
            eValues_1storder(i) = 1 / nValues(i);
            eValues_2ndorder(i) = 1 / nValues(i)^2;
            results(i, 4) = eValues_1storder(i);
            results(i, 5) = eValues_2ndorder(i);
        end
    end


    subplot(3,1,iCase)
    loglog(nValues, eValues, sty{1}), hold on;
    loglog(nValues, eValues_1storder, sty{2});
    loglog(nValues, eValues_2ndorder, sty{3});
    title(sprintf('%s-%s', BCTypeInner, BCTypeOuter));
    xlim([ nMin nMax ]);
    legend(leg);

end

disp(sprintf('#N\tDirichlet-Dirichlet\tNeumann-Dirichlet\tDirichlet-Neumann\t1storder\t2ndorder'));
for i = 1: nPoints
    disp(sprintf('%d\t%e\t%e\t%e\t%e\t%e', nValues(i), results(i,1), results(i,2), results(i,3), results(i,4), results(i,5)));
end
