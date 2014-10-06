%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
domain = [0,lengthX;[-1,1]*2.25*lengthY];
resolutionX = 500;
timespan = [0,4*lengthX/u];

%% Velocity definition
perturbationCase = 3; %@\label{ll:Bickley jet velocity definition start}
phiTimespan = [0,25];
phiInitial = [0,0];
phiSol = ode45(@d_phi,phiTimespan,phiInitial);
timeResolution = 1e5;
phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),1);
phi1Max = max(phi1);
lDerivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase,phiSol,phi1Max); %@\label{ll:Bickley jet velocity definition end}
incompressible = true;
periodicBc = [true,false]; %@\label{ll:Bickley jet periodicBc}

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);

%% LCS parameters
% Lambda-lines
lambdaRange = .8:.01:1.1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);

% Shrinklines
shrinklineMaxLength = 1e8;
shrinklineLocalMaxDistance = 8*gridSpace;
shrinklineOdeSolverOptions = odeset('relTol',1e-4);

% Stretchlines
stretchlineMaxLength = 1e8;
stretchlineLocalMaxDistance = 4*gridSpace;
stretchlineOdeSolverOptions = odeset('relTol',1e-4);

%% Lambda-line LCSs
% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [3.25e6,1.5e6;1.4e6,2.6e6]; %@\label{ll:Bickley jet Poincare start}
poincareSection(2).endPosition = [6.5e6,-1.4e6; 5.e6,-3.e6];
poincareSection(3).endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection(4).endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection(5).endPosition = [1.65e7,1.5e6;1.5e7,2.6e6]; %@\label{ll:Bickley jet Poincare end}

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(80);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

...
for lambda = lambdaRange    
...    
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);          
    closedLambdaLineCandidate = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions);        
...    
end
...   
%% Hyperbolic shrinkline LCSs
[shrinklineLcs,shrinklineLcsInitialPosition] = seed_curves_from_lambda_max(shrinklineLocalMaxDistance,shrinklineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'periodicBc',periodicBc,'odeSolverOptions',shrinklineOdeSolverOptions);

%% Hyperbolic stretchline LCSs
[stretchlineLcs,stretchlineLcsInitialPosition] = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'periodicBc',periodicBc,'odeSolverOptions',stretchlineOdeSolverOptions);
