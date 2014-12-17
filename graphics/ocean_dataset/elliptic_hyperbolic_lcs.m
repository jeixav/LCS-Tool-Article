function elliptic_hyperbolic_lcs

%% Input parameters
timespan = [100,130];
domain = [0,6;-34,-28];
resolutionX = 400;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

%% Velocity definition
load('ocean_geostrophic_velocity.mat')
% Set velocity to zero at boundaries
vLon(:,[1,end],:) = 0;
vLon(:,:,[1,end]) = 0;
vLat(:,[1,end],:) = 0;
vLat(:,:,[1,end]) = 0;
interpMethod = 'spline';
vLonInterpolant = griddedInterpolant({time,lat,lon},vLon,interpMethod);
vLatInterpolant = griddedInterpolant({time,lat,lon},vLat,interpMethod);
lDerivative = @(t,x,~)derivative(t,x,vLonInterpolant,vLatInterpolant);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgEigenvalueFromMainGrid = false;
cgAuxGridRelDelta = .01;

% Lambda lines
poincareSection = struct('endPosition',{},'numPoints',{},...
    'orbitMaxLength',{});
poincareSection(1).endPosition = [3.3,-32.1;3.7,-31.6];
poincareSection(2).endPosition = [1.3,-30.9;2.0,-31.2];
poincareSection(3).endPosition = [4.9,-29.6;5.7,-29.6];
poincareSection(4).endPosition = [4.9,-31.4;5.3,-31.4];
poincareSection(5).endPosition = [3.0,-29.3;3.5,-29.3];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),...
        diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 4*(2*pi*rOrbit);
end
lambda = .9:.02:1.1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-5,'initialStep',1e-2);
forceEtaComplexNaN = true;

% Shrink lines
shrinkLineMaxLength = 20;
shrinkLineLocalMaxDistance = 4*deltaX;
shrinkLineOdeSolverOptions = odeset('relTol',1e-4);

% Stretch lines
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 8*deltaX;
stretchLineOdeSolverOptions = odeset('relTol',1e-3);

% Graphic properties
repellingColor = 'r';
attractingColor = 'b';
ellipticColor = [0,.6,0];

hAxes = setup_figure(domain);
xlabel(hAxes,'Longitude (°)')
ylabel(hAxes,'Latitude (°)')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,...
    resolution,timespan,...
    'incompressible',incompressible,...
    'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,...
    'auxGridRelDelta',cgAuxGridRelDelta);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))

%% Elliptic LCSs
% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),...
    input.endPosition(:,2)),poincareSection,'UniformOutput',false);
hPoincareSection = [hPoincareSection{:}];
set(hPoincareSection,'color',ellipticColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',ellipticColor)
set(hPoincareSection,'MarkerEdgeColor','w')

[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(...
    domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,...
    'forceEtaComplexNaN',forceEtaComplexNaN,...
    'odeSolverOptions',lambdaLineOdeSolverOptions);

ellipticLcs = elliptic_lcs(closedLambdaLinePos);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)];

% Plot elliptic LCSs
hEllipticLcs = plot_elliptic_lcs(hAxes,ellipticLcs);
set(hEllipticLcs,'color',ellipticColor)
set(hEllipticLcs,'linewidth',2)

%% Hyperbolic repelling LCSs
shrinkLine = seed_curves_from_lambda_max(...
    shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),...
    cgEigenvector(:,1:2),domain,resolution,...
    'odeSolverOptions',shrinkLineOdeSolverOptions);

% Remove shrink lines inside elliptic LCSs
for i = 1:nPoincareSection
    shrinkLine = remove_strain_in_elliptic(shrinkLine,ellipticLcs{i});
end

% Plot hyperbolic repelling LCSs
hRepellingLcs = cellfun(@(position)plot(hAxes,position(:,1),...
    position(:,2)),shrinkLine,'UniformOutput',false);
hRepellingLcs = [hRepellingLcs{:}];
set(hRepellingLcs,'color',repellingColor)

uistack(hEllipticLcs,'top')
uistack(hPoincareSection,'top')

filename = 'elliptic_repelling_lcs.tikz';
hFigure = get(hAxes,'parent');
matlab2tikz(filename,'showInfo',false,'width','\figurewidth',...
    'figurehandle',hFigure)

%% Hyperbolic attracting LCSs
hAxes = setup_figure(domain);
xlabel(hAxes,'Longitude (°)')
ylabel(hAxes,'Latitude (°)')

% Plot finite-time Lyapunov exponent
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))

% Copy objects from repelling LCS plot
hPoincareSection = copyobj(hPoincareSection,hAxes);
hEllipticLcs = copyobj(hEllipticLcs,hAxes);

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minima
stretchLine = seed_curves_from_lambda_max(...
    stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),...
    cgEigenvector(:,3:4),domain,resolution,...
    'odeSolverOptions',stretchLineOdeSolverOptions);

% Remove stretch lines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchLine = remove_strain_in_elliptic(stretchLine,ellipticLcs{i});
end

% Plot hyperbolic attracting LCSs
hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),...
    position(:,2)),stretchLine,'UniformOutput',false);
hAttractingLcs = [hAttractingLcs{:}];
set(hAttractingLcs,'color',attractingColor)

uistack(hEllipticLcs,'top')
uistack(hPoincareSection,'top')

filename = 'elliptic_attracting_lcs.tikz';
hFigure = get(hAxes,'parent');
matlab2tikz(filename,'showInfo',false,'width','\figurewidth',...
    'figurehandle',hFigure)
