% ROUNDSHOW a 3D round show of a plot.
% (scheme package)
%
% roundshow (pause)
% -----------------
%
% a 3D round show, simply changes the view in regular intervals.
%
% Input
% -----
% - pause::value: inverse speed {DEFAULT: no pause, fast == 0}
%
% Output
% ------
% none
%
% Example
% -------
% plot_tree (sample_tree); shine ('-p -a');
% roundshow;
%
% See also
% Uses
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function roundshow (speed)

if (nargin < 1)||isempty(speed),
    speed = 0; % {DEFAULT: very fast, no pause}
end

for ward = 0 : 5 : 360,
    view([ward-37.5 30]);
    axis vis3d
    if speed ~= 0,
        pause (speed);
    end
    drawnow;
%     saveas(gcf,sprintf('roundshow%d.png',ward),'png')%'-dpng')
    
    tprint(sprintf('roundshow%0.3d.png',ward),'-HR-png',[],'-a')
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
end