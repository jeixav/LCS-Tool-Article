function hyperbolic_shear_lcs_details_pdf

lcs_tool_root = fullfile('..','..','..','..','lcs_toolbox');
demo = fullfile(lcs_tool_root,'demo','double_gyre');

oldFolder = cd(demo);

try
    addpath(fullfile('..','..'))
    hyperbolic_shear_lcs_details
    cd(oldFolder)
    
    hFigure = 1;
    filename = 'hyperbolic_shear_lcs_details_strainline';
    savefig(hFigure,filename)
    delete_title(hFigure)
    filenamePDF = [filename,'.pdf'];
    if exist(filenamePDF,'file')
        delete(filenamePDF)
    end
    print_pdf(hFigure,filename)
    
    % FIXME hFigure should not be hard-coded; should probably search for title
    % string
    hFigure = 6;
    filename = 'hyperbolic_shear_lcs_details_stretchline';
    savefig(hFigure,filename)
    delete_title(hFigure)
    filenamePDF = [filename,'.pdf'];
    if exist(filenamePDF,'file')
        delete(filenamePDF)
    end
    print_pdf(hFigure,filename)
catch err
    cd(oldFolder)
    rethrow(err)
end

function delete_title(hFigure)

hAxes = findobj(hFigure,'type','axes','-not','Tag','Colorbar');
hTitle = get(hAxes,'title');
delete(hTitle)
