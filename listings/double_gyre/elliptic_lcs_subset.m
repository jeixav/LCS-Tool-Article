%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolution = [500,250];
timespan = [0,10];

%% Velocity definition
lDerivative = @(t,x)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
cgStrainOdeSolverOptions = odeset('relTol',1e-5); %@\label{ll:double gyre cgStrainOdeSolverOptions}

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
lambdaRange = .93:.01:1.07; %@\label{ll:double gyre lambdaRange}

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

%% Lambda-line LCSs

for lambda = lambdaRange %@\label{ll:double gyre lambda loop start}
...        
    closedLambdaLineCandidate = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',showGraph);
...
end %@\label{ll:double gyre lambda loop end}
