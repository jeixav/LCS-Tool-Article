%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
timespan = [0,4*lengthX/u];
domain = [0,lengthX;[-1,1]*2.25*lengthY];
resolutionX = 500;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

%% Velocity definition
perturbationCase = 3; %@\label{ll:Bickley jet velocity definition start}
phiTimespan = [0,25];
phiInitial = [0,0];
phiSol = ode45(@d_phi,phiTimespan,phiInitial);
timeResolution = 1e5;
phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution), 1);
phi1Max = max(phi1);
lDerivative = @(t,x)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase,phiSol,phi1Max); %@\label{ll:Bickley jet velocity definition end}
incompressible = true;
periodicBc = [true,false]; %@\label{ll:Bickley jet periodicBc}

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);

%% LCS parameters
% Lambda lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [3.25e6,1.5e6;1.4e6,2.6e6]; %@\label{ll:Bickley jet Poincare start}
poincareSection(2).endPosition = [6.5e6,-1.4e6;5e6,-3e6];
poincareSection(3).endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection(4).endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection(5).endPosition = [1.65e7,1.5e6;1.5e7,2.6e6]; %@\label{ll:Bickley jet Poincare end}
[poincareSection.numPoints] = deal(20);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)), diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end
lambda = .9:.01:1.1;
lambdaLineOdeSolverOptions = odeset('absTol',10);
forceEtaComplexNaN = true;

% Shrink lines
shrinkLineMaxLength = 1e8;
shrinkLineLocalMaxDistance = 8*deltaX;
shrinkLineOdeSolverOptions = odeset('absTol',1);

% Stretch lines
stretchLineMaxLength = 1e8;
stretchLineLocalMaxDistance = 8*deltaX;
stretchLineOdeSolverOptions = odeset('absTol',10);

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);

%% Elliptic LCSs
[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,'forceEtaComplexNaN',forceEtaComplexNaN,'odeSolverOptions',lambdaLineOdeSolverOptions,'periodicBc',periodicBc);

ellipticLcs = elliptic_lcs(closedLambdaLinePos);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)];

%% Hyperbolic repelling LCSs
shrinkLine = seed_curves_from_lambda_max(shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinkLineOdeSolverOptions,'periodicBc',periodicBc);

% Remove shrink lines inside elliptic LCSs
for i = 1:nPoincareSection
    shrinkLine = remove_strain_in_elliptic(shrinkLine,ellipticLcs{i});
end

%% Hyperbolic attracting LCSs
stretchLine = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions,'periodicBc',periodicBc);

% Remove stretch lines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchLine = remove_strain_in_elliptic(stretchLine,ellipticLcs{i});
end
