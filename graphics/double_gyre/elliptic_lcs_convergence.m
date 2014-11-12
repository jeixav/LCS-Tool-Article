function elliptic_lcs_convergence

% Parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolutionX = 500:250:1000;
timespan = [0,10];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

% Lambda-lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [.55,.55;.1,.1];
poincareSection(2).endPosition = [1.53,.45;1.95,.05];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end
lambda = .93:.01:1.07;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
forceEtaComplexNaN = true;

% Graphics properties
ellipticColor = [0,.6,0];

for m = 1:numel(resolutionX)
    % Make x and y grid spacing as equal as possible
    lResolutionX = resolutionX(m);
    resolutionY = equal_resolution(domain,resolutionX);
    resolution = [lResolutionX,resolutionY];
    
    hAxes = setup_figure(domain);
    hFigure = get(hAxes,'parent');
    delete(get(hAxes,'xlabel'))
    delete(get(hAxes,'ylabel'))
    
    % Cauchy-Green strain eigenvalues and eigenvectors
    warnState = warning('off','eig_cgStrain:unequalDelta');
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);
    warning(warnState)
    
    %% Elliptic LCSs
    warnState = warning('off','integrate_line:isDiscontinuousLargeAngle');
    warning('off','integrate_line:initialConditionNaN')
    warning('off','poincare_closed_orbit:selectNeighborOrbit')
    warning('off','poincare_closed_orbit:bisectionOpenOrbit')
    [closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,'forceEtaComplexNaN',forceEtaComplexNaN,'lambdaLineOdeSolverOptions',lambdaLineOdeSolverOptions);
    warning(warnState)

    ellipticLcs = elliptic_lcs(closedLambdaLinePos);
    ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)]; %#ok<AGROW>
    
    % Plot elliptic LCSs
    hEllipticLcs = plot_elliptic_lcs(hAxes,ellipticLcs);
    set(hEllipticLcs,'color',ellipticColor)
    set(hEllipticLcs,'linewidth',2)

    % Plot closed lambda lines
    hClosedLambdaLinePos = plot_closed_orbit(hAxes,closedLambdaLinePos);
    hClosedLambdaLineNeg = plot_closed_orbit(hAxes,closedLambdaLineNeg);
    hClosedLambdaLine = [hClosedLambdaLinePos,hClosedLambdaLineNeg];
    set(hClosedLambdaLine,'color',ellipticColor)
    drawnow

    filename = strcat('elliptic_lcs_convergence_',num2str(lResolutionX),'.tikz');
    matlab2tikz(filename,'showInfo',false,'relativeDataPath',fullfile('graphics','double_gyre'),'width','\figurewidth','figurehandle',hFigure)
end
