% Strainlines LCSs
strainlineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLocalMaxDistance = 2*gridSpace; %@\label{ll:strainlineLocalMaxDistance}

strainlineLcs = seed_curves_from_lambda_max(strainlineLocalMaxDistance,strainlineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution);

% Stretchlines LCSs
stretchlineMaxLength = 20;
stretchlineLocalMaxDistance = 10*gridSpace;  %@\label{ll:stretchlineLocalMaxDistance}

stretchlineLcs = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution);
