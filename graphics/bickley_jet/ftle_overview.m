function ftle_overview

%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
domain = [0,lengthX;[-1,1]*2.25*lengthY];

%% Velocity definition
perturbationCase = 3;
phiTimespan = [0,25];
phiInitial = [0,0];
phiSol = ode45(@d_phi,phiTimespan,phiInitial);
timeResolution = 1e5;
phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),1);
phi1Max = max(phi1);
lDerivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase,phiSol,phi1Max);
incompressible = true;

% Make x and y grid spacing as equal as possible
resolutionX = 500;
gridSpace = diff(domain(1,:))/(double(resolutionX)-1);
resolutionY = round(diff(domain(2,:))/gridSpace);
resolution = [resolutionX,resolutionY];

timespan = [0,4*lengthX/u];

hAxes = setup_figure(domain);

%% Cauchy-Green strain eigenvalues and eigenvectors
cgEigenvalue = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
hFigure = 1;
delete(findobj(hFigure,'Type','axes','Tag','Colorbar'))

filename = 'ftle_overview';
print_pdf(hFigure,filename)
