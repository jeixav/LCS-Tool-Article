% Plot FTLE at various resolutions

function ftle_resolutions

lcs_tool_root = fullfile('..','..','..','..','lcs_toolbox');
demo = fullfile(lcs_tool_root,'demo','bickley_jet');

oldFolder = cd(demo);
addpath(fullfile('..','..'))

resolutionX = [100,250,500,1000];

%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
domain = [0,lengthX;[-1,1]*2.25*lengthY];
timespan = [0,4*lengthX/u];

perturbationCase = 3;
lDerivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase);
incompressible = true;

for m = 1:numel(resolutionX)
   gridSpace = diff(domain(1,:))/(double(resolutionX(m))-1);
   resolutionY = round(diff(domain(2,:))/gridSpace);
   resolution = [resolutionX(m),resolutionY];
   
   cgEigenvalue = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);
   cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
   ftle_ = ftle(cgEigenvalue2,diff(timespan));

   hAxes = setup_figure(domain);
   title(hAxes,['Resolution: (',num2str(resolution(1)),',',num2str(resolution(2)),')'])
   plot_ftle(hAxes,domain,resolution,ftle_);
   colormap(hAxes,flipud(gray))
   filename = fullfile(oldFolder,['ftle_',num2str(resolution(1)),'_',num2str(resolution(2))]);
   print_pdf(gcf,filename)
   close(gcf)
end

cd(oldFolder)
