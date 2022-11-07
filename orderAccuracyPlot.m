close all

N = [6, 11, 21, 41, 81, 161];
dx = 1./N;
errors = [6.362e-02, 6.113e-02, 1.260e-02, 2.744e-03, 5.766e-04, 1.246e-04];
slope2 = 1./N.^2;
slope3 = 1./N.^3;

figure(1);
s2p = loglog(N,10*slope2,'b--', 'LineWidth', 1.5);
hold on
s3p = loglog(N,60*slope3,'g--', 'LineWidth', 1.5);
erp = loglog(N,errors,'ro-', 'LineWidth', 2);
xlabel('N');
ylabel('Error');
legend([erp,s2p,s3p], 'Error', '2nd order', '3rd order');
