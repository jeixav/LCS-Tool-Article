% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
domain = [2e6,.55*lengthX;[-1,.25]*2.25*lengthY];
resolutionX = 500;
timespan = [0,2*lengthX/u];

% Velocity definition
perturbationCase = 3;
phiTimespan = [0,25];
phiInitial = [0,0];
phiSol = ode45(@d_phi,phiTimespan,phiInitial);
timeResolution = 1e5;
phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),1);
phi1Max = max(phi1);
lDerivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase,phiSol,phi1Max);
incompressible = true;

% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-4);

% Strainlines
strainlineMaxLength = 1e8;
strainlineOdeSolverOptions = odeset('relTol',1e-4);

% Make x and y grid spacing as equal as possible
gridSpace = diff(domain(1,:))/(double(resolutionX)-1);
resolutionY = round(diff(domain(2,:))/gridSpace) + 1;
resolution = [resolutionX,resolutionY];

strainlineLocalMaxDistance = 4*gridSpace;

% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

% Strainline LCSs
strainlineLcs = seed_curves_from_lambda_max(strainlineLocalMaxDistance,strainlineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',strainlineOdeSolverOptions);
