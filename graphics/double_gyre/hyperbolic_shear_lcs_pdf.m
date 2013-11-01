function hyperbolic_shear_lcs_pdf

lcs_tool_root = fullfile('..','..','..','..','lcs_toolbox');
double_gyre_demo = fullfile(lcs_tool_root,'demo','double_gyre');

addpath(lcs_tool_root)
addpath(double_gyre_demo)

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
