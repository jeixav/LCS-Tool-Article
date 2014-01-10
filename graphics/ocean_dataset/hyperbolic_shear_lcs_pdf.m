function hyperbolic_shear_lcs_pdf

hyperbolic_shear_lcs

hFigure = 1;
filename = 'hyperbolic_shear_lcs_strainline';
delete_title(hFigure)
print_pdf(hFigure,filename)

hFigure = 2;
filename = 'hyperbolic_shear_lcs_stretchline';
delete_title(hFigure)
print_pdf(hFigure,filename)

function delete_title(hFigure)

hAxes = get(hFigure,'children');
hTitle = get(hAxes,'title');
delete(hTitle)
