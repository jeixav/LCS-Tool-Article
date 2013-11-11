function hyperbolic_shear_lcs_pdf

lcs_tool_root = fullfile('..','..','..','..','lcs_toolbox');
demo = fullfile(lcs_tool_root,'demo','ocean_dataset');

oldFolder = cd(demo);

try
    addpath(fullfile('..','..'))
    hyperbolic_shear_lcs
    cd(oldFolder)
    
    hFigure = 1;
    filename = 'hyperbolic_shear_lcs_forward';
    savefig(hFigure,filename)
    delete_title(hFigure)
    filenamePDF = [filename,'.pdf'];
    if exist(filenamePDF,'file')
        delete(filenamePDF)
    end
    print_pdf(hFigure,filename)
    
    hFigure = 2;
    filename = 'hyperbolic_shear_lcs_backward';
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

hAxes = get(hFigure,'children');
hTitle = get(hAxes,'title');
delete(hTitle)
