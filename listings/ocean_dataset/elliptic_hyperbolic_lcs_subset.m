% Input parameters
domain = [0,6;-34,-28];
res = [400,400]; %@\label{ll:ocean data resolution}
time = [100,130]; %@\label{ll:ocean data timespan}

% Velocity definition
load('ocean_geostrophic_velocity.mat');
...
intMet = 'spline'; %@\label{ll:interpolation start}
vLonI = griddedInterpolant({time,lat,lon},vLon,intMet);
vLatI = griddedInterpolant({time,lat,lon},vLat,intMet);
lDerivative = @(t,x,~)flowdata_derivative(t,x,vLonI,vLatI); %@\label{ll:interpolation end}

% LCS parameters
% Cauchy-Green strain
cgVmg = false; %@\label{ll:cgEigenvalueFromMainGrid}
cgAgd = .01; %@\label{ll:cgAuxGridRelDelta}

% Lambda lines
...
ps(1).endPosition = [3.3,-32.1;3.7,-31.6]; %@\label{ll:ocean data poincareSection(1)}
...
lambda = .9:.02:1.1; %@\label{ll:ocean data lambdaRange}
llOptions = odeset('relTol',1e-6);
...
% Shrink lines
shrinkLMaxLength = 20;
shrinkLMaxDistance = 2*gridSpace;
...
% Stretch lines
stretchLMaxLength = 20;
stretchLMaxDistance = 4*gridSpace;
...  
%% Cauchy-Green strain eigenvalues and eigenvectors
[cgV,cgD] = eig_cgStrain(lDerivative,domain,res,time, 'eigenvalueFromMainGrid',cgVmg, 'auxGridRelDelta',cgAgd); %@\label{ll:eig_cgStrain}
%% Elliptic LCSs
[closedLlp,closedLln] = poincare_closed_orbit_range(domain,res,cgV,cgD, lambda,ps,'odeSolverOptions',llOptions);%@\label{ll:ocean data poincare_closed_orbit_multi}
...
% Hyperbolic shrink line LCSs
shrinkL = seed_curves_from_lambda_max( shrinkLMaxDistance,shrinkLMaxLength,cgV(:,2), cgD(:,1:2),domain,res); %@\label{ll:ocean data shrinkLineLcs}
...
% Hyperbolic stretch line LCSs
stretchL = seed_curves_from_lambda_max( stretchLMaxDistance,stretchLMaxLength,-cgV(:,1), cgD(:,3:4),domain,res); %@\label{ll:ocean data stretchLineLcs}
