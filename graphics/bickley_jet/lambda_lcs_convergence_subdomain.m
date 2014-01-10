function lambda_lcs_convergence_subdomain

%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
domain = [2e6,lengthX/2;[-1,.25]*2.25*lengthY];
resolutionX = 400:100:600;
timespan = [0,2*lengthX/u];

%% Velocity definition
perturbationCase = 3;
phiTimespan = [0,25];
phiInitial = [0,0];
phiSol = ode45(@d_phi,phiTimespan,phiInitial);
timeResolution = 1e5;
phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),1);
phi1Max = max(phi1);
lDerivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase,phiSol,phi1Max);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-4);

% Lambda-lines
lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
poincareSection.endPosition = [6.5,-1.4;4.5,-3.5]*1e6;
[poincareSection.numPoints] = deal(100);
rOrbit = hypot(diff(poincareSection.endPosition(:,1)),diff(poincareSection.endPosition(:,2)));
poincareSection.orbitMaxLength = 2*(2*pi*rOrbit);
lambdaLineColor = [0,.6,0];
dThresh = 1e-3;

warning('off','eig_cgStrain:unequalDelta')
warning('off','eig_cgStrain:unequalAuxGridDelta')
warning('off','integrate_line:isDiscontinuousLargeAngle')

for m = 1:numel(resolutionX)
    % Make x and y grid spacing as equal as possible
    lResolutionX = resolutionX(m);
    gridSpace = diff(domain(1,:))/(double(lResolutionX)-1);
    resolutionY = round(diff(domain(2,:))/gridSpace) + 1;
    resolution = [lResolutionX,resolutionY];
    disp(['Resolution: ',num2str(resolution)])
    
    hAxes = setup_figure(domain);
    title(hAxes,['Resolution: ',num2str(resolution(1)),'\times',num2str(resolution(2))])
    
    %% Cauchy-Green strain eigenvalues and eigenvectors
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);
    
    % Plot finite-time Lyapunov exponent
    cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
    ftle_ = ftle(cgEigenvalue2,diff(timespan));
    plot_ftle(hAxes,domain,resolution,ftle_);
    colormap(hAxes,flipud(gray))
    drawnow

    %% Lambda-line LCSs
    % Plot Poincare sections
    hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
    set(hPoincareSection,'color',lambdaLineColor)
    set(hPoincareSection,'LineStyle','--')
    set(hPoincareSection,'marker','o')
    set(hPoincareSection,'MarkerFaceColor',lambdaLineColor)
    set(hPoincareSection,'MarkerEdgeColor','k')
    drawnow
    
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
    [closedLambdaLine,~,hFigure] = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'dThresh',dThresh,'showGraph',true);
    delete(hFigure(2))
    hPoincare = findobj(hFigure(1),'type','axes','Tag',[]);
    title(hPoincare,['Poincare return map, resolution: ',num2str(resolution(1)),'\times',num2str(resolution(2))])
    set(hPoincare,'ylim',[-2,2]*1e4)
    hLegend = findobj(hFigure(1),'type','axes','Tag','legend');
    set(hLegend,'Location','SouthWest')
    
    % Plot lambda-line LCSs
    hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
    hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
    hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
    set(hLambdaLineLcs,'color',lambdaLineColor)
    set(hLambdaLineLcs,'linewidth',2)
    
    % Plot all closed lambda lines
    hClosedLambdaLinePos = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{1}{1});
    hClosedLambdaLineNeg = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{1}{2});
    hClosedLambdaLine = horzcat(hClosedLambdaLinePos,hClosedLambdaLineNeg);
    set(hClosedLambdaLine,'color',lambdaLineColor)
    uistack(hPoincareSection,'top')
    drawnow
end
