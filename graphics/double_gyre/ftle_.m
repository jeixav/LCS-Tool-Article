function ftle_

%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolutionX = 500;
timespan = [0,5];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

% Make x and y grid spacing as equal as possible
gridSpace = diff(domain(1,:))/(double(resolutionX)-1);
resolutionY = round(diff(domain(2,:))/gridSpace);
resolution = [resolutionX,resolutionY];
    
hAxes = setup_figure(domain);
title(hAxes,['Resolution: ',num2str(resolution(1)),'\times',num2str(resolution(2))])
    
%% Cauchy-Green strain eigenvalues and eigenvectors
s = warning('off','eig_cgStrain:unequalDelta');
cgEigenvalue = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);
warning(s)

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
delete(get(hAxes,'title'))
hFigure = get(hAxes,'parent');
print_pdf(hFigure,'ftle.pdf')
