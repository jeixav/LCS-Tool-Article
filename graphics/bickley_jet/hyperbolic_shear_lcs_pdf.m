function hyperbolic_shear_lcs_pdf

lcs_tool_root = fullfile('..','..','..','..','lcs_toolbox');
bickley_jet_demo = fullfile(lcs_tool_root,'demo','bickley_jet');

addpath(lcs_tool_root)
addpath(bickley_jet_demo)

hyperbolic_shear_lcs

hFigure = 1;
filename = 'hyperbolic_shear_lcs_forward';
savefig(hFigure,filename)

delete_title(hFigure)

if exist(filename,'file')
    delete(filename)
end
print_pdf(hFigure,filename)

hFigure = 2;
filename = 'hyperbolic_shear_lcs_backward';
savefig(hFigure,filename)

delete_title(hFigure)

if exist(filename,'file')
    delete(filename)
end
print_pdf(hFigure,filename)

function delete_title(hFigure)

hAxes = get(hFigure,'children');
hTitle = get(hAxes,'title');
delete(hTitle)
