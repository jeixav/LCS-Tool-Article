%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
domain = [2e6,lengthX/2;[-1,.25]*2.25*lengthY];
resolutionX = 400:100:600;
timespan = [0,2*lengthX/u];

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

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-4);

% Lambda-lines
lambda = .995;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
poincareSection.endPosition = [6.5,-1.4;4.5,-3.5]*1e6;
[poincareSection.numPoints] = deal(100);
rOrbit = hypot(diff(poincareSection.endPosition(:,1)),diff(poincareSection.endPosition(:,2)));
poincareSection.orbitMaxLength = 2*(2*pi*rOrbit);
lambdaLineColor = [0,.6,0];
dThresh = 1e-4;

for m = 1:numel(resolutionX)
    % Make x and y grid spacing as equal as possible
    lResolutionX = resolutionX(m);
    gridSpace = diff(domain(1,:))/(double(lResolutionX)-1);
    resolutionY = round(diff(domain(2,:))/gridSpace) + 1;
    resolution = [lResolutionX,resolutionY];
   
    %% Cauchy-Green strain eigenvalues and eigenvectors
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);
    
    %% Lambda-line LCSs
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
    [closedLambdaLine,~,hFigure] = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'dThresh',dThresh,'showGraph',true);
end
