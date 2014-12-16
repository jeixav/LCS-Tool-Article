commit d1700275fbb35a7c16481b37cb6b61e53810e421
Author: Kristjan Onu <jeixav@gmail.com>
Date:   Tue Dec 16 17:00:00 2014 -0500

    Update Bickley jet demo parameters

diff --git a/graphics/bickley_jet/elliptic_hyperbolic_lcs.m b/graphics/bickley_jet/elliptic_hyperbolic_lcs.m
index 7329b70..141ae28 100644
--- a/graphics/bickley_jet/elliptic_hyperbolic_lcs.m
+++ b/graphics/bickley_jet/elliptic_hyperbolic_lcs.m
@@ -5,8 +5,8 @@ u = 62.66;
 lengthX = pi*earthRadius;
 lengthY = 1.77e6;
 epsilon = [.075,.4,.3];
-domain = [0,lengthX;[-1,1]*2.25*lengthY];
 timespan = [0,4*lengthX/u];
+domain = [0,lengthX;[-1,1]*2.25*lengthY];
 resolutionX = 500;
 [resolutionY,deltaX] = equal_resolution(domain,resolutionX);
 resolution = [resolutionX,resolutionY];
@@ -17,51 +17,56 @@ phiTimespan = [0,25];
 phiInitial = [0,0];
 phiSol = ode45(@d_phi,phiTimespan,phiInitial);
 timeResolution = 1e5;
-phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),1);
+phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),...
+    timeResolution),1);
 phi1Max = max(phi1);
-lDerivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase,phiSol,phi1Max);
+lDerivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,...
+    perturbationCase,phiSol,phi1Max);
 incompressible = true;
 periodicBc = [true,false];
 
 %% LCS parameters
-
-% Lambda-lines
-poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
+% Lambda lines
+poincareSection = struct('endPosition',{},'numPoints',{},...
+    'orbitMaxLength',{});
 poincareSection(1).endPosition = [3.25e6,1.5e6;1.4e6,2.6e6];
 poincareSection(2).endPosition = [6.5e6,-1.4e6;5e6,-3e6];
 poincareSection(3).endPosition = [1e7,1.5e6;8e6,2.6e6];
 poincareSection(4).endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
 poincareSection(5).endPosition = [1.65e7, 1.5e6;1.5e7, 2.6e6];
-[poincareSection.numPoints] = deal(80);
+[poincareSection.numPoints] = deal(20);
 nPoincareSection = numel(poincareSection);
 for i = 1:nPoincareSection
-    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
+    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),...
+        diff(poincareSection(i).endPosition(:,2)));
     poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
 end
 lambda = .9:.01:1.1;
-lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
+lambdaLineOdeSolverOptions = odeset('absTol',10);
 forceEtaComplexNaN = true;
 
-% Strainlines
-strainlineMaxLength = 1e8;
-strainlineLocalMaxDistance = 8*deltaX;
-strainlineOdeSolverOptions = odeset('relTol',1e-4);
+% Shrink lines
+shrinkLineMaxLength = 1e8;
+shrinkLineLocalMaxDistance = 8*deltaX;
+shrinkLineOdeSolverOptions = odeset('absTol',1);
 
-% Stretchlines
-stretchlineMaxLength = 1e8;
-stretchlineLocalMaxDistance = 4*deltaX;
-stretchlineOdeSolverOptions = odeset('relTol',1e-4);
+% Stretch lines
+stretchLineMaxLength = 1e8;
+stretchLineLocalMaxDistance = 8*deltaX;
+stretchLineOdeSolverOptions = odeset('absTol',10);
 
-% Graphics properties
+% Graphic properties
 repellingColor = 'r';
 attractingColor = 'b';
 ellipticColor = [0,.6,0];
-lcsInitialPositionMarkerSize = 2;
+initialPositionMarkerSize = 2;
 
 hAxes = setup_figure(domain);
 
 %% Cauchy-Green strain eigenvalues and eigenvectors
-[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);
+[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,...
+    resolution,timespan,...
+    'incompressible',incompressible);
 
 % Plot finite-time Lyapunov exponent
 cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
@@ -71,7 +76,8 @@ colormap(hAxes,flipud(gray))
 
 %% Elliptic LCSs
 % Plot Poincare sections
-hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection,'UniformOutput',false);
+hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),...
+    input.endPosition(:,2)),poincareSection,'UniformOutput',false);
 hPoincareSection = [hPoincareSection{:}];
 set(hPoincareSection,'color',ellipticColor)
 set(hPoincareSection,'LineStyle','--')
@@ -79,7 +85,11 @@ set(hPoincareSection,'marker','o')
 set(hPoincareSection,'MarkerFaceColor',ellipticColor)
 set(hPoincareSection,'MarkerEdgeColor','w')
 
-[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,'forceEtaComplexNaN',forceEtaComplexNaN,'odeSolverOptions',lambdaLineOdeSolverOptions,'periodicBc',periodicBc);
+[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(...
+    domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,...
+    'forceEtaComplexNaN',forceEtaComplexNaN,...
+    'odeSolverOptions',lambdaLineOdeSolverOptions,...
+    'periodicBc',periodicBc);
 
 ellipticLcs = elliptic_lcs(closedLambdaLinePos);
 ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)];
@@ -96,25 +106,34 @@ hClosedLambdaLine = [hClosedLambdaLinePos,hClosedLambdaLineNeg];
 set(hClosedLambdaLine,'color',ellipticColor)
 
 %% Hyperbolic repelling LCSs
-[strainlineLcs,strainlineLcsInitialPosition] = seed_curves_from_lambda_max(strainlineLocalMaxDistance,strainlineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',strainlineOdeSolverOptions,'periodicBc',periodicBc);
+[shrinkLine,shrinkLineInitialPosition] = seed_curves_from_lambda_max(...
+    shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),...
+    cgEigenvector(:,1:2),domain,resolution,...
+    'odeSolverOptions',shrinkLineOdeSolverOptions,...
+    'periodicBc',periodicBc);
 
-% Remove strainlines inside elliptic LCSs
+% Remove shrink lines inside elliptic LCSs
 for i = 1:nPoincareSection
-    strainlineLcs = remove_strain_in_elliptic(strainlineLcs,ellipticLcs{i});
-    idx = inpolygon(strainlineLcsInitialPosition(1,:),strainlineLcsInitialPosition(2,:),ellipticLcs{i}(:,1),ellipticLcs{i}(:,2));
-    strainlineLcsInitialPosition = strainlineLcsInitialPosition(:,~idx);
+    shrinkLine = remove_strain_in_elliptic(shrinkLine,ellipticLcs{i});
+    idx = inpolygon(shrinkLineInitialPosition(1,:),...
+        shrinkLineInitialPosition(2,:),ellipticLcs{i}(:,1),...
+        ellipticLcs{i}(:,2));
+    shrinkLineInitialPosition = shrinkLineInitialPosition(:,~idx);
 end
 
 % Plot hyperbolic repelling LCSs
-hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs,'UniformOutput',false);
-hStrainlineLcs = [hStrainlineLcs{:}];
-set(hStrainlineLcs,'color',repellingColor)
-hStrainlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineLcsInitialPosition(1,idx),strainlineLcsInitialPosition(2,idx)),1:size(strainlineLcsInitialPosition,2),'UniformOutput',false);
-hStrainlineLcsInitialPosition = [hStrainlineLcsInitialPosition{:}];
-set(hStrainlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
-set(hStrainlineLcsInitialPosition,'marker','o')
-set(hStrainlineLcsInitialPosition,'MarkerEdgeColor','w')
-set(hStrainlineLcsInitialPosition,'MarkerFaceColor',repellingColor)
+hRepellingLcs = cellfun(@(position)plot(hAxes,position(:,1),...
+    position(:,2)),shrinkLine,'UniformOutput',false);
+hRepellingLcs = [hRepellingLcs{:}];
+set(hRepellingLcs,'color',repellingColor)
+hShrinkLineInitialPosition = arrayfun(@(idx)plot(hAxes,...
+    shrinkLineInitialPosition(1,idx),shrinkLineInitialPosition(2,idx)),...
+    1:size(shrinkLineInitialPosition,2),'UniformOutput',false);
+hShrinkLineInitialPosition = [hShrinkLineInitialPosition{:}];
+set(hShrinkLineInitialPosition,'MarkerSize',initialPositionMarkerSize)
+set(hShrinkLineInitialPosition,'marker','o')
+set(hShrinkLineInitialPosition,'MarkerEdgeColor','w')
+set(hShrinkLineInitialPosition,'MarkerFaceColor',repellingColor)
 
 uistack(hEllipticLcs,'top')
 uistack(hClosedLambdaLine,'top')
@@ -122,7 +141,8 @@ uistack(hPoincareSection,'top')
 
 filename = 'elliptic_repelling_lcs.tex';
 hFigure = get(hAxes,'parent');
-matlab2tikz(filename,'showInfo',false,'width','\figurewidth','figurehandle',hFigure)
+matlab2tikz(filename,'showInfo',false,'width','\figurewidth',...
+    'figurehandle',hFigure)
 
 %% Hyperbolic attracting LCSs
 hAxes = setup_figure(domain);
@@ -138,31 +158,36 @@ hEllipticLcs = copyobj(hEllipticLcs,hAxes);
 
 % FIXME Part of calculations in seed_curves_from_lambda_max are
 % unsuitable/unecessary for stretchlines do not follow ridges of λ₁
-% minimums
-[stretchlineLcs,stretchlineLcsInitialPosition] = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineOdeSolverOptions,'periodicBc',periodicBc);
-
-% FIXME Discard stretchlines to prevent pdflatex from running out of memory
-idx = 1:3:numel(stretchlineLcs);
-stretchlineLcs = stretchlineLcs(idx);
-stretchlineLcsInitialPosition = stretchlineLcsInitialPosition(:,idx);
-
-% Remove stretchlines inside elliptic LCSs
+% minima
+[stretchLine,stretchLineInitialPosition] = seed_curves_from_lambda_max(...
+    stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),...
+    cgEigenvector(:,3:4),domain,resolution,...
+    'odeSolverOptions',stretchLineOdeSolverOptions,...
+    'periodicBc',periodicBc);
+
+% Remove stretch lines inside elliptic LCSs
 for i = 1:nPoincareSection
-    stretchlineLcs = remove_strain_in_elliptic(stretchlineLcs,ellipticLcs{i});
-    idx = inpolygon(stretchlineLcsInitialPosition(1,:),stretchlineLcsInitialPosition(2,:),ellipticLcs{i}(:,1),ellipticLcs{i}(:,2));
-    stretchlineLcsInitialPosition = stretchlineLcsInitialPosition(:,~idx);
+    stretchLine = remove_strain_in_elliptic(stretchLine,ellipticLcs{i});
+    idx = inpolygon(stretchLineInitialPosition(1,:),...
+        stretchLineInitialPosition(2,:),ellipticLcs{i}(:,1),...
+        ellipticLcs{i}(:,2));
+    stretchLineInitialPosition = stretchLineInitialPosition(:,~idx);
 end
 
 % Plot hyperbolic attracting LCSs
-hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs,'UniformOutput',false);
-hStretchlineLcs = [hStretchlineLcs{:}];
-set(hStretchlineLcs,'color',attractingColor)
-hStretchlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineLcsInitialPosition(1,idx),stretchlineLcsInitialPosition(2,idx)),1:size(stretchlineLcsInitialPosition,2),'UniformOutput',false);
-hStretchlineLcsInitialPosition = [hStretchlineLcsInitialPosition{:}];
-set(hStretchlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
-set(hStretchlineLcsInitialPosition,'marker','o')
-set(hStretchlineLcsInitialPosition,'MarkerEdgeColor','w')
-set(hStretchlineLcsInitialPosition,'MarkerFaceColor',attractingColor)
+hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),...
+    position(:,2)),stretchLine,'UniformOutput',false);
+hAttractingLcs = [hAttractingLcs{:}];
+set(hAttractingLcs,'color',attractingColor)
+hStretchLineInitialPosition = arrayfun(@(idx)plot(hAxes,...
+    stretchLineInitialPosition(1,idx),...
+    stretchLineInitialPosition(2,idx)),...
+    1:size(stretchLineInitialPosition,2),'UniformOutput',false);
+hStretchLineInitialPosition = [hStretchLineInitialPosition{:}];
+set(hStretchLineInitialPosition,'MarkerSize',initialPositionMarkerSize)
+set(hStretchLineInitialPosition,'marker','o')
+set(hStretchLineInitialPosition,'MarkerEdgeColor','w')
+set(hStretchLineInitialPosition,'MarkerFaceColor',attractingColor)
 
 uistack(hEllipticLcs,'top')
 uistack(hClosedLambdaLine,'top')
@@ -170,4 +195,5 @@ uistack(hPoincareSection,'top')
 
 filename = 'elliptic_attracting_lcs.tex';
 hFigure = get(hAxes,'parent');
-matlab2tikz(filename,'showInfo',false,'width','\figurewidth','figurehandle',hFigure)
+matlab2tikz(filename,'showInfo',false,'width','\figurewidth',...
+    'figurehandle',hFigure)
