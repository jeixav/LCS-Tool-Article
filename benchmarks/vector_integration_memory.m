% Use: "profile -memory on" to measure memory use
% Reference: http://undocumentedmatlab.com/blog/undocumented-profiler-options/
%
% profile('-memory','on')
% vector_integration_memory
% stats = profile('info');
% strcmp('vector_integration_memory',stats.FunctionTable(1).FunctionName)
% stats.FunctionTable(1).TotalTime

function vector_integration_memory

% Parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolutionX = 1000;
timespan = [0,10];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

resolutionY = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];
    
[~,~] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);
