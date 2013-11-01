function hyperbolic_shear_lcs_pdf

lcs_tool_root = fullfile('..','..','..','..','lcs_toolbox');
ocean_dataset_demo = fullfile(lcs_tool_root,'demo','ocean_dataset');

addpath(lcs_tool_root)
addpath(ocean_dataset_demo)

hyperbolic_shear_lcs

hFigure = 1;
filename = 'hyperbolic_shear_lcs_forward';
savefig(hFigure,filename)

delete_title(hFigure)

if exist([filename,'.pdf'],'file')
    delete([filename,'.pdf'])
end
print_pdf(hFigure,filename)

hFigure = 2;
filename = 'hyperbolic_shear_lcs_backward';
savefig(hFigure,filename)

delete_title(hFigure)

if exist([filename,'.pdf'],'file')
    delete([filename,'.pdf'])
end
print_pdf(hFigure,filename)

function delete_title(hFigure)

hAxes = get(hFigure,'children');
hTitle = get(hAxes,'title');
delete(hTitle)
