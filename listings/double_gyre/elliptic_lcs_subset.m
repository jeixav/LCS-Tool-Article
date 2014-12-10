%% Input parameters
...
domain = [0,2;0,1];
res = [500,250];
time = [0,10];

%% Velocity definition
lDerivative = @(t,x)derivative(t,x,$\epsilon$,A,$\omega$);

%% LCS parameters
cgOptions = odeset('relTol',1e-5); %@\label{ll:double gyre cgStrainOdeSolverOptions}

% Lambda lines
...
ps(1).endPosition = [.55,.55;.1,.1]; %@\label{ll:double gyre poincareSection(1)}
ps(2).endPosition = [1.53,.45;1.95,.05]; %@\label{ll:double gyre poincareSection(2)}
...
lambda = .93:.01:1.07; %@\label{ll:double gyre lambdaRange}
llOptions = odeset('relTol',1e-6); %@\label{ll:double gyre lambdaLineOdeSolverOptions}

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgV,cgD] = eig_cgStrain(lDerivative,domain,res,time, 'odeSolverOptions',cgOptions);

%% Elliptic LCSs
[closedLlp,closedLln] = poincare_closed_orbit_range(domain,res,cgV,cgD, lambda,ps,'odeSolverOptions',llOptions); %@\label{ll:double gyre lambda range}
ellipticLcs = elliptic_lcs(closedLlp);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLln)]; %@\label{ll:double gyre elliptic_lcs}
