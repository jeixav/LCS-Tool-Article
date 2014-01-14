% print_pdf Save figure to PDF
%
% Syntax
% print_pdf(figure,filename)
%
% Description
% The output of saveas(figure,'figure.pdf','pdf') is ugly! This function
% sets a few graphics properties to improve the appearance and then uses
% epstopdf to produce a PDF file
%
% Example
% x = -pi:.1:pi;
% y = sin(x);
% plot(x,y)
% print_pdf(gcf,'sin.pdf')

function print_pdf(figure,filename)

narginchk(2,2)

axes = findobj(figure,'type','axes');

if numel(axes) > 1
    hTitle = cell2mat(get(axes,'title'));
    hXlabel = cell2mat(get(axes,'xlabel'));
    hYlabel = cell2mat(get(axes,'ylabel'));
else
    hTitle = get(axes,'title');
    hXlabel = get(axes,'xlabel');
    hYlabel = get(axes,'ylabel');
end

allHandles = [axes;hTitle;hXlabel;hYlabel];

originalFontname = get(allHandles,'fontname');
originalFontsize = get(allHandles,'fontsize');
originalPaperPositionMode = get(figure,'PaperPositionMode');

pdfFontname = 'lucida';
pdfFontsize = 14;
pdfPaperPositionMode = 'auto';

set(allHandles,'fontname',pdfFontname)
set(allHandles,'fontsize',pdfFontsize)
set(figure,'PaperPositionMode',pdfPaperPositionMode)

[pathstr,name,ext] = fileparts(filename);

if isempty(ext)
    ext = 'pdf';
else
    % Remove . from start of ext
    ext = ext(2:end);
end

epsFullfile = fullfile(tempdir,[name,'.eps']);
if exist(epsFullfile,'file')
    error(['EPS file exists: ',epsFullfile])
end
saveas(figure,epsFullfile,'epsc2')

pdfFullfile = fullfile(pathstr,[name,'.',ext]);
if exist(pdfFullfile,'file')
    delete(epsFullfile)
    error(['PDF file exists: ',pdfFullfile])
end

[status,~] = unix(['convert -density 300 ',epsFullfile,' ',pdfFullfile]);
if status
    error([mfilename,':unixConvert'],'Unix convert command error')
end
delete(epsFullfile)

allHandlesC = mat2cell(allHandles,ones(length(allHandles),1));
cellfun(@(handle,value)set(handle,'fontname',value),allHandlesC,originalFontname)
cellfun(@(handle,value)set(handle,'fontsize',value),allHandlesC,originalFontsize)
set(figure,'PaperPositionMode',originalPaperPositionMode)
