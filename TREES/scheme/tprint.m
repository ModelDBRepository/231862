% TPRINT   Simplified printing.
% (scheme package)
%
% [name path] = tprint (name, options, fsize, options)
% ----------------------------------------------------
%
% prints the current figure to a file.
%
% Input
% -----
% - name::string: name of the output-file without extension
%     {DEFAULT : open gui fileselect} spaces and other weird symbols not
%     allowed!
% - options::string: {DEFAULT '-R -jpg'}
%    '-SHR'  1200 dpi
%    '-HR'   600  dpi
%    '-R'    300  dpi
%    '-LR'   150  dpi
%    '-jpg'  as jpg
%    '-tif'  as tif
%    '-eps'  as eps  (more than one output format is possible!)
% - fsize::2-tupel: fixed size in cm [horiz. x vertical] {DEFAULT 15cm x 10cm}
% - options::string: {DEFAULT: ''} , see shine
%    '-a'  : axis invisible
%    '-f'  : full axis
%    '-s'  : scalebar in um
%    '-p'  : camlight and lighting gouraud
%    '-3d' : 3D view
%
% Output
% ------
% - name::string: name of output file; [] no file was selected -> no output
% - path::sting: path of the file, complete string is therefore: [path name]
%
% Example
% -------
% plot_tree (sample_tree); tprint;
%
% See also
% Uses shine
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function [tname, tpath] = tprint (tname, options, fsize, shineoptions)

if (nargin < 1)||isempty(tname)
    [tname, tpath] = uiputfile ('.*','print figure content','fig');
    if tname == 0
        tname = [];
        return
    end
else
    tpath = '';
end
nstart = unique ([0 strfind(tname, '/') strfind(tname, '\')]);
if nstart (end) > 0
    tpath = [tpath tname(1 : nstart (end))];
    tname(1 : nstart (end)) = '';
end

tname(strfind(tname,'%')) = 'p'; % percent symbol makes problems with ghostscript

if (nargin<2)||isempty(options)
    options = '-R';
    switch tname(end-2:end)
        case {'pdf','eps','png','tif','svg'}
            options = strcat(options,' -',tname(end-2:end));
        otherwise
            options = strcat(options,' -jpg');
    end
end

if (nargin<3)||isempty(fsize)
    %     fsize = [15 10];
    %     fsize = [21 29.7];
    set(gcf,'Units','centimeter')
    fsize = get(gcf,'Position');
    fsize = fsize(3:4);
    %     set(gcf,'Units','pixels')
end

if (nargin<4)||isempty(shineoptions)
    shineoptions = 'none';
end

if     strfind (options, '-HR')
    res = 600;
elseif strfind (options, '-SHR')
    res = 1200;
elseif strfind (options, '-LR')
    res = 150;
else
    res = 300;
end

set (gcf, 'units', 'centimeters', 'PaperUnits', 'centimeters', ...
    'papersize', fsize, 'PaperPosition', [0 0 fsize]);
set(gcf,'units','pixels')

if ~strcmp (shineoptions, 'none')
    shine (shineoptions);
end
usedFormats = 0;
if ~isempty(strfind (options, '-png'))
    if ~isempty(strfind(tname,'.png'))
        print ('-dpng', ['-r' num2str(res)], [tpath tname]);
    else
        print ('-dpng', ['-r' num2str(res)], [tpath tname '.png']);
    end
    usedFormats = usedFormats +1;
end
if ~isempty(strfind (options, '-jpg'))
    print ('-djpeg', ['-r' num2str(res)], [tpath tname '.jpg']);
    usedFormats = usedFormats +1;
end
if ~isempty(strfind(options,'-tif'))
    if ~isempty(strfind(tname,'.tif'))
        print ('-dtiff', ['-r' num2str(res)], [tpath tname]);
    else
        print ('-dtiff', ['-r' num2str(res)], [tpath tname '.tif']);
    end
    usedFormats = usedFormats +1;
end
if ~isempty(strfind(options,'-eps'))
    %     print ('-depsc2','-painters', [tpath tname '.eps']);
    print ('-painters','-depsc2', [tpath tname '.eps']);
    usedFormats = usedFormats +1;
end
if ~isempty(strfind(options,'-svg'))
    %     hold all
    %     set(gcf,'Units','data')
    %     get(gcf,'Position')
    %     line(
    % line
    axes('Position',[.005 .005 .99 .99],'xtick',[],'ytick',[],'box','on');
    set(gcf,'children',flipud(get(gcf,'children')))
    plot2svg([tpath tname '.svg'])
    %     print ('-dsvg','-painters', [tpath tname '.svg']);
    usedFormats = usedFormats +1;
end
if ~isempty(strfind(options,'-pdf'))
    if ~isempty(strfind(tname,'.pdf'))
        print ('-painters','-dpdf', ['-r' num2str(res)],'-loose', [tpath tname]);
    else
        print ('-painters','-dpdf', ['-r' num2str(res)], [tpath tname '.pdf']);
    end
    usedFormats = usedFormats +1;
end
if ~usedFormats
    error('No file format definition found in option string "%s"',options)
end