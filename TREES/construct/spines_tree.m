% SPINES_TREE   Add spines to an existing tree.
% (trees package)
%
% tree = spines_tree (intree, XYZ, dneck, dhead, mlneck, stdlneck, ipart,
% options)
% -------------------------------------------------------------------------
%
% attaches cylinders with diameter dhead (dhead = length) to closest node
% on an existing tree, introducing a neck with diameter dneck. If XYZ
% coordinates are not defined the spine is attached to a randomly picked
% node with distance mlneck+-stdlneck. If region with name "spines" exists
% then nodes are appended to that region otherwise region named "spines" is
% created.
%
% Input
% -----
% - intree::integer:index of tree in trees or structured tree
% - XYZ:: matrix [X Y Z] or just a number of spines to add {DEFAULT: 100}
% - dneck::value: diameter of spine neck {DEFAULT: 0.25 um}
% - dhead::value: diameter of spine head {DEFAULT: 1 um}
% - mlneck::value: mean neck length {DEFAULT: 1 um}
% - stdlneck::value: standard deviation neck length {DEFAULT: 1 um}
% - ipart::index:index to the subpart to be considered {DEFAULT: all nodes}
%     (needs to be real index values not logical subset)
% - options::string: {DEFAULT '-w'}
%     '-s' : show
%     '-w' : waitbar
%     '-sr': instead of just adding a region "spines" it adds region for
%            spine heads and necks
%
% Output
% ------
% if no output is declared the trees are added in trees
% - tree:: structured output tree
%
% Example
% -------
% spines_tree (sample_tree, 300, 0.25, 1, 1, 1, [], '-s')
%
% See also quaddiameter_tree MST_tree
% Uses 
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function [tree, indhead, indneck] = spines_tree (intree, XYZ, dneck, dhead, mlneck, stdlneck, ipart, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length (trees); % {DEFAULT tree: last tree in trees cell array}
end;

ver_tree (intree); % verify that input is a tree structure

% use full tree for this function
if ~isstruct (intree),
    tree = trees {intree};
else
    tree = intree;
end
if nargin < 8 || isempty(options)
    options = '';
end
X = tree.X; % X-locations of nodes on tree
Y = tree.Y; % Y-locations of nodes on tree
Z = tree.Z; % Z-locations of nodes on tree

N = size (X, 1); % number of nodes in tree

if (nargin<5)||isempty(mlneck),
    mlneck = 1;
end

if (nargin<6)||isempty(stdlneck),
    stdlneck = 1;
end

if (nargin<2)||isempty(XYZ),
    XYZ = 100;
end

if (nargin<3)||isempty(dneck),
    dneck = .5;
end

if (nargin<4)||isempty(dhead),
    dhead = 1;
end

if (nargin < 7)||isempty(ipart),
    ipart = (1 : N)'; % {DEFAULT index: select all nodes/points}
end


    dir = dir_tree(tree,'-n');
    if numel(XYZ)==1,
        indy = ceil (rand (XYZ, 1) * length (ipart));
    elseif all(XYZ<N) % they are indices
        indy = XYZ;
    end
    dXYZ = zeros(numel(indy),3);
    for n = 1:numel(indy)
        [~,~,V] = svd(dir(ipart(indy(n)),:),0);  % calculate an orthogonal vector to the direction of the dendrite
        M = makehgtform('axisrotate',dir(ipart(indy(n)),:),2*rand()*pi); %Rotate the vector around axis dendrite axis by random radians
        dXYZ(n,:) = V(:,3)' * M(1:3,1:3);       
    end
    XYZ  = [X(ipart (indy)) Y(ipart (indy)) Z(ipart (indy))] + ...
        repmat((randn (numel(indy), 1) * stdlneck + mlneck),1,3) .* dXYZ;

      
    
if isfield (tree, 'R')
    if isfield (tree, 'rnames'),
        
        if ~isempty(strfind(options,'-sr'))
            r = find (strcmpi (tree.rnames, 'spine_neck'));
            if ~isempty (r);
                iR(1) = r (1);
            else
                iR(1) = max( max(tree.R),numel(tree.rnames)) + 1;
                flag = 1;
            end
            r = find (strcmpi (tree.rnames, 'spine_head'));
            if ~isempty (r);
                iR(2) = r (1);
            else
                iR(2) = max( max(tree.R),numel(tree.rnames)) + 1 + flag;
            end
        else
            iR = find (strcmp (tree.rnames, 'spines'));
            if ~isempty (iR);
                iR = iR (1);
            else
                iR = max(tree.R,numel(tree.rnames)) + 1;
                tree.rnames = [tree.rnames, 'spines'];
            end
        end
    else
        iR = max(tree.R,numel(tree.rnames)) + 1;
    end
else
    iR = 1;
end

if numel(iR) == 1
    iR(2) = iR(1);
end
if ~isempty(strfind(options,'-sr'))
    tree.rnames{iR(1)} = 'spines_neck';
    tree.rnames{iR(2)} = 'spines_head';
end



if (nargin<8)||isempty(options),
    options = '-w';
end
if strfind (options, '-w'), % waitbar option: initialization
    HW = waitbar (0, 'spining...');
    set (HW, 'Name', '..PLEASE..WAIT..YEAH..');
end
for ward = 1 : size (XYZ, 1),
    if strfind (options, '-w'), % waitbar option: update
        if mod (ward, 100) == 0,
            waitbar (ward / size (XYZ, 1), HW);
        end
    end
    [tree,indneck] = insert_tree (tree, ...
        [1 iR(1) XYZ(ward, 1) XYZ(ward, 2) XYZ(ward, 3) dneck ipart(indy(ward))], 'none');
    [tree,indhead] = insert_tree (tree, ...
        [1 iR(2) XYZ(ward, 1)+dXYZ(ward,1)*dhead, ...
        XYZ(ward, 2)+dXYZ(ward,2)*dhead XYZ(ward, 3)+dXYZ(ward,3)*dhead dhead N+1+2*(ward - 1)], 'none');
end
if strfind (options, '-w'), % waitbar option: close
    close (HW);
end

if strfind (options, '-s'),
    clf; shine; hold on;
    plot_tree (tree);
    title  ('add spines to your tree');
    xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
    view (3); grid on; axis image;
end

if (nargout == 1)||(isstruct(intree))
    varargout {1}  = tree; % if output is defined then it becomes the tree
else
    trees {intree} = tree; % otherwise the orginal tree in trees is replaced
end

