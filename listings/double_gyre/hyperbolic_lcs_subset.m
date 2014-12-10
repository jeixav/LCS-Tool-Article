% Shrink lines
shrinkLMaxLength = 20;
shrinkLMaxDistance = 2*gridSpace; %@\label{ll:shrinkLineLocalMaxDistance}

% Stretch lines
stretchLMaxLength = 20;
stretchLMaxDistance = 10*gridSpace; %@\label{ll:stretchLineLocalMaxDistance}

%% Hyperbolic shrink line LCSs
shrinkL = seed_curves_from_lambda_max( shrinkLMaxDistance,shrinkLMaxLength,cgV(:,2), cgD(:,1:2),domain,res);
...
%% Hyperbolic stretch line LCSs
stretchL = seed_curves_from_lambda_max( stretchLMaxDistance,stretchLMaxLength,-cgV(:,1), cgD(:,3:4),domain,res);
...
