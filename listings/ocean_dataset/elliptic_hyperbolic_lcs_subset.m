% Input parameters
domain = [0,6;-34,-28];
resolution = [400,400]; %@\label{ll:ocean resolution}
timespan = [100,130]; %@\label{ll:ocean timespan}

% Velocity definition
load('ocean_geostrophic_velocity.mat');
...
interpMethod = 'spline'; %@\label{ll:interpolation start}
vLonInterpolant = griddedInterpolant({time,lat,lon},vLon,interpMethod);
vLatInterpolant = griddedInterpolant({time,lat,lon},vLat,interpMethod);
lDerivative = @(t,x)flowdata_derivative(t,x,vLonInterpolant,vLatInterpolant); %@\label{ll:interpolation end}
incompressible = true; %@\label{ll:incompressible}

% LCS parameters
% Cauchy-Green strain
cgEigenvalueFromMainGrid = false; %@\label{ll:cgEigenvalueFromMainGrid}
cgAuxGridRelDelta = .01; %@\label{ll:cgAuxGridRelDelta}

% Lambda lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [3.3,-32.1;3.7,-31.6]; %@\label{ll:ocean data poincareSection(1)}
poincareSection(2).endPosition = [1.3,-30.9;2.0,-31.2];
poincareSection(3).endPosition = [4.9,-29.6;5.7,-29.6];
poincareSection(4).endPosition = [4.9,-31.4;5.3,-31.4];
poincareSection(5).endPosition = [3.0,-29.3;3.5,-29.3]; %@\label{ll:ocean data poincareSection(5)}
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 4*(2*pi*rOrbit);
end
lambda = .9:.02:1.1; %@\label{ll:ocean lambda}
lambdaLineOdeSolverOptions = odeset('relTol',1e-6,'initialStep',1e-2);
forceEtaComplexNaN = true;

% Shrink lines
shrinkLineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
shrinkLineLocalMaxDistance = 2*gridSpace;
shrinkLineOdeSolverOptions = odeset('relTol',1e-4);

% Stretch lines
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 4*gridSpace;
stretchLineOdeSolverOptions = odeset('relTol',1e-4);
...
% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta); %@\label{ll:eig_cgStrain}

%% Elliptic LCSs
[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,'forceEtaComplexNaN',forceEtaComplexNaN,'odeSolverOptions',lambdaLineOdeSolverOptions); %@\label{ll:ocean elliptic LCS start}

ellipticLcs = elliptic_lcs(closedLambdaLinePos);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)]; %@\label{ll:ocean elliptic LCS end}

% Hyperbolic shrink line LCSs
shrinkLine = seed_curves_from_lambda_max(shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinkLineOdeSolverOptions); %@\label{ll:ocean hyperbolic start}

% Remove shrink lines inside elliptic LCSs
for i = 1:nPoincareSection
    shrinkLine = remove_strain_in_elliptic(shrinkLine,ellipticLcs{i});
end

% Hyperbolic stretch line LCSs
stretchLine = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

% Remove stretch lines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchLine = remove_strain_in_elliptic(stretchLine,ellipticLcs{i});
end %@\label{ll:ocean hyperbolic end}
