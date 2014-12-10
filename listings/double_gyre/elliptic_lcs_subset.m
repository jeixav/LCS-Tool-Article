%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolutionX = 500;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];
timespan = [0,10];

%% Velocity definition
lDerivative = @(t,x)derivative(t,x,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
cgOptions = odeset('relTol',1e-5); %@\label{ll:double gyre cgStrainOdeSolverOptions}

% Lambda-lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [.55,.55;.1,.1]; %@\label{ll:double gyre poincareSection(1)}
poincareSection(2).endPosition = [1.53,.45;1.95,.05]; %@\label{ll:double gyre poincareSection(2)}
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end
lambdaLineOdeSolverOptions = odeset('relTol',1e-6); %@\label{ll:double gyre lambdaLineOdeSolverOptions}
lambda = .93:.01:1.07; %@\label{ll:double gyre lambdaRange}

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

%% Elliptic LCSs
[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,ps,'odeSolverOptions',lambdaLineOdeSolverOptions); %@\label{ll:double gyre lambda range}
ellipticLcs = elliptic_lcs(closedLambdaLinePos);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)]; %@\label{ll:double gyre elliptic_lcs}
