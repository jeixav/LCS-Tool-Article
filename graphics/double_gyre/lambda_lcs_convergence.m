function lambda_lcs_convergence

%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolutionX = 500:250:1000;
timespan = [0,5];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

% Lambda-lines
lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [.55,.55;.2,.5];
poincareSection(2).endPosition = [1.53,.45;1.9,.5];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(100);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

% Graphics properties
lambdaLineColor = [0,.6,0];

for m = 1:numel(resolutionX)
    % Make x and y grid spacing as equal as possible
    lResolutionX = resolutionX(m);
    gridSpace = diff(domain(1,:))/(double(lResolutionX)-1);
    resolutionY = round(diff(domain(2,:))/gridSpace);
    resolution = [lResolutionX,resolutionY];
    disp(['Resolution: ',num2str(resolution)])
    
    hAxes = setup_figure(domain);
    title(hAxes,['Resolution: ',num2str(resolution(1)),'\times',num2str(resolution(2))])
    
    %% Cauchy-Green strain eigenvalues and eigenvectors
    s = warning('off','eig_cgStrain:unequalDelta');
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);
    warning(s)
    
    %% Lambda-line LCSs
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
    warnState = warning('off','integrate_line:isDiscontinuousLargeAngle');
    closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions);
    warning(warnState)
    
    % Plot lambda-line LCSs
    hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
    hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
    hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
    set(hLambdaLineLcs,'color',lambdaLineColor)
    set(hLambdaLineLcs,'linewidth',2)
    
    % Plot all closed lambda lines
    hClosedLambdaLinePos = cell(nPoincareSection,1);
    hClosedLambdaLineNeg = cell(nPoincareSection,1);
    for j = 1:nPoincareSection
        hClosedLambdaLinePos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{1});
        hClosedLambdaLineNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{2});
    end
    hClosedLambdaLine = horzcat(hClosedLambdaLinePos{:},hClosedLambdaLineNeg{:});
    set(hClosedLambdaLine,'color',lambdaLineColor)
    drawnow
    hFigure = get(hAxes,'parent');
    print_pdf(hFigure,['lambda_lcs_convergence_',num2str(resolution(1))])
end
