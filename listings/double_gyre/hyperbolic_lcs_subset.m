% Shrink lines
shrinkLineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
shrinkLineLocalMaxDistance = 2*gridSpace; %@\label{ll:shrinkLineLocalMaxDistance}
shrinkLineOdeSolverOptions = odeset('relTol',1e-6);

% Stretch lines
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 10*gridSpace; %@\label{ll:stretchLineLocalMaxDistance}
stretchLineOdeSolverOptions = odeset('relTol',1e-6);

%% Hyperbolic shrink line LCSs
[shrinkLineLcs,shrinkLineLcsInitialPosition] = seed_curves_from_lambda_max(shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinkLineOdeSolverOptions);
...
%% Hyperbolic stretch line LCSs
[stretchLineLcs,stretchLineLcsInitialPosition] = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);
...