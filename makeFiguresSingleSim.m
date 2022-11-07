% Make figures that present the results in different ways

close all

x_vector = linspace(0, L_x, NI);
y_vector = linspace(0, L_y, NJ);
z_vector = linspace(0, L_z, NK);
theta = linspace(0, 2*pi, 1000);
x_outerCircle = outerRadius * cos(theta) + xCenter;
y_outerCircle = outerRadius * sin(theta) + yCenter;
x_innerCircle = innerRadius * cos(theta) + xCenter;
y_innerCircle = innerRadius * sin(theta) + yCenter;
kIndex = ceil(NK/2);
nIsoLines = 20;

figure(1);
bitmapPlot = imagesc(x_vector, y_vector, T(:,:,kIndex)');
axis equal;
set(gca,'YDir','normal');
set(bitmapPlot, 'AlphaData', ~isnan(T(:,:,kIndex)'));
colorbar;
title('Bitmap of Temperature');
xlabel('x');
ylabel('y');
hold on;
plot(x_innerCircle, y_innerCircle, 'r-', 'LineWidth', 1.5);
plot(x_outerCircle, y_outerCircle, 'r-', 'LineWidth', 1.5);

figure(2);
[colors, contourPlot] = contourf(x_vector, y_vector, T(:,:,kIndex)', nIsoLines);
axis equal;
%set(contourPlot, 'LineColor', 'none');
colorbar;
title('Contour plot of Temperature');
xlabel('x');
ylabel('y');
hold on;
plot(x_innerCircle, y_innerCircle, 'r-', 'LineWidth', 2);
plot(x_outerCircle, y_outerCircle, 'r-', 'LineWidth', 2);

figure(3);
surface(x_vector, y_vector, T(:,:,kIndex)');
axis equal;
axis([x_vector(1), x_vector(end), y_vector(1), y_vector(end)]);
colorbar;
shading interp;
title('Surface plot of Temperature');
xlabel('x');
ylabel('y');
hold on;
%plot(x_innerCircle, y_innerCircle, 'r-', 'LineWidth', 2);
%plot(x_outerCircle, y_outerCircle, 'r-', 'LineWidth', 2);

figure(4);
errorAbs(ghostNodesIndices) = NaN;
[colors, contourPlot] = contourf(x_vector, y_vector, errorAbs(:,:,kIndex)', nIsoLines);
axis equal;
%set(contourPlot, 'LineColor', 'none');
colorbar;
title('Contour plot of error distribution');
xlabel('x');
ylabel('y');
hold on;
plot(x_innerCircle, y_innerCircle, 'r-', 'LineWidth', 2);
plot(x_outerCircle, y_outerCircle, 'r-', 'LineWidth', 2);

figure(5);
surface(x_vector, y_vector, errorAbs(:,:,kIndex)');
axis equal;
axis([x_vector(1), x_vector(end), y_vector(1), y_vector(end)]);
colorbar;
shading interp;
title('Surface plot of error distribution');
xlabel('x');
ylabel('y');
hold on;
%plot(x_innerCircle, y_innerCircle, 'r-', 'LineWidth', 2);
%plot(x_outerCircle, y_outerCircle, 'r-', 'LineWidth', 2);

if caseType == "Cyl"
    zPlane = L_z/2;
elseif caseType == "Sphere"
    zPlane = zCenter;
end

figure(6);
slice(x_vector, y_vector, z_vector, permute(T,[2,1,3]), xCenter, yCenter, zPlane);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
colorbar;
title('Slice plot of Temperature');

figure(7);
slice(x_vector, y_vector, z_vector, permute(errorAbs,[2,1,3]), xCenter, yCenter, zPlane);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
colorbar;
title('Slice plot of error');

figure(8);
semilogy(normHistory/normHistory(1), 'LineWidth', 2);
title("Development of the 2-norm of the change between time levels (convergence)");
ylabel("$\displaystyle\frac{||T^n-T^{n-1}||}{||T^1-T^0||}$",'interpreter','latex');
xlabel("n, time level");

if (caseType == "Cyl") % If cylinder case, compare with Martin's code
    
    addpath("./ehsan gamle filer/");
    qNeumann = 0;
    if bcType == "DirNeu"
        bcTypeOuter = 'Neumann';
        qNeumann = tempGradOuter;
    else
        bcTypeOuter = 'Dirichlet';
    end
    if bcType == "NeuDir"
        bcTypeInner = 'Neumann';
        qNeumann = tempGradInner;
    else
        bcTypeInner = 'Dirichlet';
    end
    [L2_err, T_implicit, error_implicit] = laplace_Martin(NI, bcTypeOuter, bcTypeInner, 0, T_inner, T_outer, qNeumann, xCenter, yCenter, innerRadius, outerRadius, false);
    
    figure(9);
    [colors, contourPlot] = contourf(x_vector, y_vector, T_implicit', nIsoLines);
    axis equal;
    %set(contourPlot, 'LineColor', 'none');
    colorbar;
    title('Contour plot of Temperature, implicit code');
    xlabel('x');
    ylabel('y');
    hold on;
    plot(x_innerCircle, y_innerCircle, 'r-', 'LineWidth', 2);
    plot(x_outerCircle, y_outerCircle, 'r-', 'LineWidth', 2);
    
    figure(10);
    [colors, contourPlot] = contourf(x_vector, y_vector, error_implicit', nIsoLines);
    axis equal;
    %set(contourPlot, 'LineColor', 'none');
    colorbar;
    title('Contour plot of error, implicit code');
    xlabel('x');
    ylabel('y');
    hold on;
    plot(x_innerCircle, y_innerCircle, 'r-', 'LineWidth', 2);
    plot(x_outerCircle, y_outerCircle, 'r-', 'LineWidth', 2);
    
    figure(11);
    [colors, contourPlot] = contourf(x_vector, y_vector, (T(:,:,kIndex)-T_implicit)', nIsoLines);
    axis equal;
    %set(contourPlot, 'LineColor', 'none');
    colorbar;
    title('Contour plot of difference between transient and implicit code');
    xlabel('x');
    ylabel('y');
    hold on;
    plot(x_innerCircle, y_innerCircle, 'r-', 'LineWidth', 2);
    plot(x_outerCircle, y_outerCircle, 'r-', 'LineWidth', 2);
    
end

spreadFigures;


