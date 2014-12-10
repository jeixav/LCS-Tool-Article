%% Input parameters
...
domain = [0,$L_x$;[-1,1]*2.25*$L_y$];
res = [500,276];
time = [0,4*$L_x$/u];

%% Velocity definition
...
lDerivative = @(t,x)derivative(t,x,u,$L_x$,$L_y$,$\epsilon$, perturbationCase,phiSol,phi1Max); %@\label{ll:Bickley jet velocity definition}
pBc = [true,false]; %@\label{ll:Bickley jet periodicBc}

%% LCS parameters
% Lambda lines
lambda = .8:.01:1.1;
...
% Shrink lines
shrinkLMaxLength = 1e8;
shrinkLMaxDistance = 8*gridSpace;
...
% Stretch lines
stretchLMaxLength = 1e8;
stretchLMaxDistance = 4*gridSpace;
...
%% Lambda line LCSs
...
ps(1).endPosition = [3.25e6,1.5e6;1.4e6,2.6e6]; %@\label{ll:Bickley jet Poincare start}
ps(2).endPosition = [6.5e6,-1.4e6; 5.e6,-3.e6];
ps(3).endPosition = [1e7,1.5e6;8e6,2.6e6];
ps(4).endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
ps(5).endPosition = [1.65e7,1.5e6;1.5e7,2.6e6]; %@\label{ll:Bickley jet Poincare end}
...
%% Cauchy-Green strain eigenvalues and eigenvectors
[cgV,cgD] = eig_cgStrain(lDerivative,domain,res,time);
%% Elliptic LCSs
[closedLlp,closedLln] = poincare_closed_orbit_range( domain,res,cgV,cgD,lambda,ps,'periodicBc',pBc);
...
%% Hyperbolic shrink line LCSs
shrinkL = seed_curves_from_lambda_max(shrinkLMaxDistance, shrinkLMaxLength,cgV(:,2),cgD(:,1:2),domain,res, 'periodicBc',pBc);
...
%% Hyperbolic stretch line LCSs
stretchL = seed_curves_from_lambda_max(stretchLMaxDistance, stretchLMaxLength,-cgV(:,1),cgD(:,3:4),domain,res, 'periodicBc',pBc);
...
