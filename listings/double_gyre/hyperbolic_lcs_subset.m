% Shrinklines
shrinklineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
shrinklineLocalMaxDistance = 2*gridSpace; %@\label{ll:shrinklineLocalMaxDistance}
strainlineOdeSolverOptions = odeset('relTol',1e-6);

% Stretchlines
stretchlineMaxLength = 20;
stretchlineLocalMaxDistance = 10*gridSpace; %@\label{ll:stretchlineLocalMaxDistance}
stretchlineOdeSolverOptions = odeset('relTol',1e-6);

%% Hyperbolic shrinkline LCSs
[shrinklineLcs,shrinklineLcsInitialPosition] = seed_curves_from_lambda_max(shrinklineLocalMaxDistance,shrinklineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinklineOdeSolverOptions);
...
%% Hyperbolic stretchline LCSs
[stretchlineLcs,stretchlineLcsInitialPosition] = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineOdeSolverOptions);
...