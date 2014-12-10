% Shrink lines
shrinkLineMaxLength = 20;
shrinkLineLocalMaxDistance = 2*deltaX; %@\label{ll:shrinkLineLocalMaxDistance}
shrinkLineOdeSolverOptions = odeset('relTol',1e-6);

% Stretch lines
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 10*deltaX; %@\label{ll:stretchLineLocalMaxDistance}
stretchLineOdeSolverOptions = odeset('relTol',1e-6);

%% Hyperbolic repelling LCSs
shrinkLine = seed_curves_from_lambda_max(shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinkLineOdeSolverOptions);

% Remove shrink lines inside elliptic LCSs
for i = 1:nPoincareSection
    shrinkLine = remove_strain_in_elliptic(shrinkLine,ellipticLcs{i});
end

%% Hyperbolic attracting LCSs
stretchLine = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

% Remove stretch lines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchLine = remove_strain_in_elliptic(stretchLine,ellipticLcs{i});
end
