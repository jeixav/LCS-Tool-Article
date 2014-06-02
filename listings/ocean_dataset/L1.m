% Input parameters
domain = [0,6;-34,-28];
resolution = [400,400];
timespan = [100,130];

% Velocity definition
load('ocean_geostrophic_velocity.mat');
...
interpMethod = 'spline';
vlon_interpolant = griddedInterpolant({time,lat,lon},vlon,interpMethod);
vlat_interpolant = griddedInterpolant({time,lat,lon},vlat,interpMethod);
lDerivative = @(t,x,~)flowdata_derivative(t,x,vlon_interpolant,vlat_interpolant);
incompressible = true;

% LCS parameters
% Cauchy-Green strain
cgEigenvalueFromMainGrid = false;
cgAuxGridRelDelta = 0.01;

% Lambda-lines
lambdaStep = 0.02; lambdaRange = 0.90:lambdaStep:1.10;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);

% Shrinklines
shrinklineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
shrinklineLocalMaxDistance = 2*gridSpace;
shrinklineOdeSolverOptions = odeset('relTol',1e-4);

% Stretchlines
stretchlineMaxLength = 20;
stretchlineLocalMaxDistance = 4*gridSpace;
stretchlineOdeSolverOptions = odeset('relTol',1e-4);
...
% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);

% Lambda-line LCSs
% Define Poincare sections; ...
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
...
poincareSection(1).endPosition = [3.3,-32.1;3.7,-31.6];
poincareSection(2).endPosition = [1.3,-30.9;1.9,-31.1];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(100);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

for lambda = lambdaRange
...
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
    closedLambdaLineCandidate = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',showGraph);
...
end
...
% Hyperbolic shrinkline LCSs
shrinklineLcs = seed_curves_from_lambda_max(shrinklineLocalMaxDistance,shrinklineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinklineOdeSolverOptions);
...
% Hyperbolic stretchline LCSs
stretchlineLcs = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineOdeSolverOptions);