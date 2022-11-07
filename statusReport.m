
% Write status report to screen, with info about progression etc.
fprintf("Sim case: %s, %s, NI=%i \n", caseType, bcType, NI);
fprintf("Simulation time: t = %.2e \n", t);
fprintf("Wall clock time: %s \n", duration(0,0,toc));
normRelative = normHistory(end)/normHistory(1);
progress = -log10(normRelative)/-log10(convCrit)*100;
fprintf("Convergence: %.2e -> %.2e  (%.1f %%) \n", normRelative, convCrit, progress);
secETA = max(0, reportInterval/(progress-previousProgress)*(100-progress));
ETA = duration(0,0,secETA);
fprintf("ETA: %s \n\n", ETA);
previousProgress = progress;
nextReportTime = nextReportTime + reportInterval;
