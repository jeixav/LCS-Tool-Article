function hyperbolic_shear_lcs_details_pdf

hyperbolic_shear_lcs_details

hFigure = 1;
filename = 'hyperbolic_shear_lcs_details_strainline';
delete_title(hFigure)
print_pdf(hFigure,filename)

% FIXME hFigure should not be hard-coded; should probably search for
% matching title string
hFigure = 12;
filename = 'hyperbolic_shear_lcs_details_stretchline';
delete_title(hFigure)
print_pdf(hFigure,filename)

function delete_title(hFigure)

hAxes = findobj(hFigure,'type','axes','-not','Tag','Colorbar');
hTitle = get(hAxes,'title');
delete(hTitle)
