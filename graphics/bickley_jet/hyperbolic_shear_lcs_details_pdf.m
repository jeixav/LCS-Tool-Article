function hyperbolic_shear_lcs_details_pdf

lcs_tool_root = fullfile('..','..','..','..','lcs_toolbox');
demo = fullfile(lcs_tool_root,'demo','bickley_jet');

addpath(lcs_tool_root,demo)

hyperbolic_shear_lcs_details
    
hFigure = 1;
delete_title(hFigure)
filename = 'hyperbolic_shear_lcs_details_strainline';
print_pdf(hFigure,filename)
    
% FIXME hFigure number should not be hard-coded; should probably search for
% title string
hFigure = 12;
delete_title(hFigure)
filename = 'hyperbolic_shear_lcs_details_stretchline';
print_pdf(hFigure,filename)

function delete_title(hFigure)

hAxes = findobj(hFigure,'type','axes','-not','Tag','Colorbar');
hTitle = get(hAxes,'title');
delete(hTitle)
