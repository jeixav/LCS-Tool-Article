function hyperbolic_shear_lcs_details_pdf

lcs_tool_root = fullfile('..','..','..','..','lcs_toolbox');
double_gyre_demo = fullfile(lcs_tool_root,'demo','double_gyre');

addpath(lcs_tool_root)
addpath(double_gyre_demo)

hyperbolic_shear_lcs_details

hFigure = 1;
delete_title(hFigure)

filename = 'hyperbolic_shear_lcs_details_forward.pdf';
if exist(filename,'file')
    delete(filename)
end
print_pdf(hFigure,filename)

% FIXME hFigure should not be hard-coded; should probably search for title
% string
hFigure = 6;
delete_title(hFigure)

filename = 'hyperbolic_shear_lcs_details_backward.pdf';
if exist(filename,'file')
    delete(filename)
end
print_pdf(hFigure,filename)

function delete_title(hFigure)

hAxes = findobj(hFigure,'type','axes','-not','Tag','Colorbar');
hTitle = get(hAxes,'title');
delete(hTitle)
