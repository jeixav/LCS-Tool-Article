function poincare_return_map

%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolution = [750,375];
timespan = [0,5];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

% Lambda-lines
poincareSection.endPosition = [.55,.55;.2,.5];
poincareSection.numPoints = 100;
rOrbit = hypot(diff(poincareSection.endPosition(:,1)),diff(poincareSection.endPosition(:,2)));
poincareSection.integrationLength = [0,2*(2*pi*rOrbit)];

lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

etaPos = lambda_line(cgEigenvector,cgEigenvalue,lambda);
s = warning('off','integrate_line:isDiscontinuousLargeAngle');
[~,~,hPoincareMap] = poincare_closed_orbit(domain,resolution,etaPos,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',true);
warning(s)

delete(findobj(hPoincareMap,'type','axes','Tag','legend'))
hAxes = findobj(hPoincareMap,'type','axes');
delete(get(hAxes,'title'))
set(hAxes,'ylim',[-1,1]*1e-4)
% Change colors to black
set(findobj(hAxes,'Marker','none'),'color',[0,0,0])
set(findobj(hAxes,'Marker','o'),'MarkerEdgeColor',[0,0,0])
set(findobj(hAxes,'Marker','o','MarkerFaceColor',[1,0,0]),'MarkerFaceColor',[0,0,0])
set(findobj(hAxes,'Marker','o','MarkerFaceColor',[0,0,1]),'MarkerFaceColor','none')

print_pdf(hPoincareMap,'poincare_return_map')
