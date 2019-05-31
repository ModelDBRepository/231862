% ROT_TREE   Rotate a tree.
% (trees package)
%
% tree = rot_tree (intree, DEG, options)
% --------------------------------------
%
% rotates a cell's anatomy by multiplying each point in space with a simple
% 3by3 (3D) or 2x2 (2D) rotation-matrix. Rotation along principal
% components is also possible. Then first pc is in x, second pc is in y,
% third pc is in z. It helps to center the tree with tran_tree beforehand
% except if rotation around a different point is required.
%
% Input
% -----
% - intree::integer:index of tree in trees or structured tree
% - DEG::1/3-tupel: degrees of rotation around resp. axis (0-360). The
%     sequence is as defined in rotation_matrix: x then y then z-axis. If
%     1-tupel rotation is in XY plane. {DEFAULT: [0 0 90]}
% - options::string: {DEFAULT: ''}
%     '-s'    : show before and after
%     '-pc2d' : directs in principal components, 2-dimensional
%     '-pc3d' : directs in principal components, 3-dimensional
%     '-m3dX' : mean axis, 3-dimensional, central axis lays on x-axis
%     '-m3dY' : mean axis, 3-dimensional, central axis lays on y-axis
%     '-m3dZ' : mean axis, 3-dimensional, central axis lays on z-axis
%     '-al'   : align region borders
%     m3d and al was implemented by Marcel Beining, 22 March 2012
%     NOTE! pc implementation does not work!
%     For the case of m3d and pc rotation first input DEG becomes index of
%     subset of nodes to be used for obtaining PCs
%
% Output
% -------
% if no output is declared the tree is changed in trees
% - tree:: structured output tree
%
% Example
% -------
% rot_tree (sample_tree, [0 30 0], '-s')
%
% See also tran_tree scale_tree flip_tree
% Uses ver_tree X Y Z
%
% -pc requires statistics toolbox for "princomp"
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function varargout = rot_tree (intree, DEG, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length(trees); % {DEFAULT tree: last tree in trees cell array}
end;

ver_tree (intree); % verify that input is a tree structure

% use full tree for this function
if ~isstruct (intree),
    tree = trees {intree};
else
    tree = intree;
end

if (nargin < 3)||isempty(options),
    options = ''; % {DEFAULT: no option}
end

if (nargin < 2)||isempty(DEG),
    if regexpi (options, '-pc'),
        DEG = (1 : length (tree.X))'; % {DEFAULT: use all nodes for PCs}
    else
        DEG = [0 0 90]; % {DEFAULT: rotate 90deg in z}
    end
end

if regexpi (options, '-pc'),
    if regexpi (options, '-pc2d'),
        XY = [tree.X tree.Y];
        p = princomp (XY (DEG, :));
        while sum (diag (p)) ~= size (p, 1),
            XY = XY * p';
            p = princomp (XY (DEG, :));
        end
        tree.X = XY (:, 1);
        tree.Y = XY (:, 2);
    else
        XYZ = [tree.X tree.Y tree.Z];
        p = princomp (XYZ (DEG, :));
        while sum (diag (p)) ~= size (p, 1),
            XYZ = XYZ * p';
            p = princomp (XYZ (DEG, :));
        end
        tree.X = XYZ (:, 1);
        tree.Y = XYZ (:, 2);
        tree.Z = XYZ (:, 3);
    end
elseif regexpi(options, '-m3d')
    % define axis to which tree is aligned: 1=x 2=y 3=z
    if     regexpi(options, '-m3dX')
        e = [0 1 0];
        raxis = 1;
        d = [2 3];
    elseif regexpi(options, '-m3dY')
%         e = [0 0 1];
        e = [1 0 0];
        raxis = 2;
        d = [3 1];
    else
        e = [1 0 0];
        raxis = 3;
        d = [1 2];
    end
    XYZ0       = [tree.X(1) tree.Y(1) tree.Z(1)];
    tree       = tran_tree (tree); % translate tree to coordinate origin
    a = zeros(1,3);
    a(raxis) = 1;
    angout = zeros(3,3);
    raxon = find(strcmpi(tree.rnames,'axon'));
    if isempty(raxon)
        raxon = 0;
    end
    for ward   = 1 : 2      % first rotate to x axis
        mXYZ   = [mean(tree.X(tree.R~=raxon)) mean(tree.Y(tree.R~=raxon)) mean(tree.Z(tree.R~=raxon))];
        mXYZ (d (1)) = 0; % make 2D norm vector orthogonal to current rotation axis
        mXYZ   = mXYZ / sqrt (sum (mXYZ.^2));
        rangle = zeros (3, 1);
        % simple angle calculation by dot product and corrections for left
        % hand rotation
        rangle(d(1)) = sign(sum(cross(a,mXYZ))) * acosd (dot (a, mXYZ));    % cross product checks if mXYZ is in clockwise direction of axis vector or not
        angout(ward,d(1)) = rangle(d(1));
        tree   = rot_tree (tree, rangle);
        d      = fliplr (d); % for next turn rotate around the 2nd axis
    end
    
%     if raxis > 1 % if tree should be rotated to axis other than x, do that
%         rangle = zeros(3,1);
%         rangle(d(raxis-1)) = sign(raxis-2.5) * 90;  % define left hand rotation around corresponding axis
%         tree = rot_tree (tree,rangle);
%     end

    % rotate tree to have maximum variation in xy (axis x and y) or xz
    % (axis z) direction using principal component analysis
    XYZ = [tree.X(tree.R~=raxon) tree.Y(tree.R~=raxon) tree.Z(tree.R~=raxon)];
    XYZ(:,raxis)= 0;    % delete information along aligned axis
    [coeff,latent]=eigs(cov(XYZ),3);
    [~,ind] = sort(diag(latent),'descend');
    coeff = coeff(:,ind);
    qual = latent(1)/latent(2) - 1; % defining quality of alignment (zero is bad because dendrites distributed equally in 2 dim)
    rangle = zeros(3,1);
%     rangle(raxis)=  sign(coeff(sigcos,1)) * acosd(dot(coeff(:,1),e));
    ang(1) = real(sign(sum(cross(e,coeff(:,1)))) * acosd(dot(coeff(:,1),e)));
    ang(2) = real(sign(sum(cross(-e,coeff(:,1)))) * acosd(dot(coeff(:,1),-e)));
% %     ang = ang + [0 180 -180];
    [~, angout(3,raxis)] = min(abs(ang));
    angout(3,raxis) = ang(angout(3,raxis));
%     rangle(raxis)=  min(acosd(dot(coeff(:,1),e)),180-acosd(dot(coeff(:,1),e))); %choose smallest angle to desired axis
%     if round(rangle(raxis))== -180
%         rangle(raxis) = 180 + rangle(raxis);
%     elseif round(rangle(raxis))== 180
%         rangle(raxis) = rangle(raxis) - 180;
%     end
    rangle(raxis) = angout(3,raxis);
    tree= rot_tree(tree,rangle);    % rotate tree to have maximum variation onto that axis
    rangle(raxis) = 0;
    
    if ~isempty(regexpi(options,'-al')) && raxis < 3 % only for y or x alignment..
        if numel(options) >= regexpi(options,'-al')+3  % get index of region to which tree should be aligned
            ind = str2num(options(regexpi(options,'-al')+3));
        else
            ind = 5;
        end
        idpar = idpar_tree(tree);   % get parent node indices
        in = find(tree.R == ind  & tree.R(idpar) == ind-1); % find points at border of region
        if isempty(in)
            warning('No points between region %s and %s found. Alignment could not be done',tree.rnames{ind},tree.rnames{ind-1})
        else
            allpoints = [tree.X(in) tree.Y(in) tree.Z(in)]; % creacte a plane of points
            allp_mean = mean(allpoints,1);
            allp_subsmean = bsxfun(@minus,allpoints,allp_mean); % Subtract "mean" point
            [~,~,V] = svd(allp_subsmean,0);
            n = V(:,3)';  % get vector orthogonal to plane
            n(setdiff(d,3)) = [];
            n = n/norm(n);
            if n(1) < 0
                n = -n;
            end
            %         figure;plot_tree(tree,tree.R)
            rangle(setdiff(d,3)) = sign(n(2)) * acosd(dot([1 0],n,2));  % get angle between plane and axis vector
            angout(:,end+1) = rangle(setdiff(d,3));  % correct tree rotation for that angle
            %     ang2(2) = acosd(dot([0 1],-n,2));
            tree = rot_tree(tree,rangle);
        end
    end
    tree       = tran_tree (tree, XYZ0);
    
else
    if numel (DEG) == 1,
        RM =  [cos(DEG) -sin(DEG);sin(DEG) cos(DEG)]; % rotation matrix 2D
        RXY = [tree.X tree.Y] * RM;
        tree.X = RXY (:, 1);
        tree.Y = RXY (:, 2);
    else
        if length (DEG) == 2,
            DEG = [DEG 0];
        end
        % rotation matrix 3D see "rotation_matrix"
        RM = rotation_matrix (deg2rad( DEG (1)), deg2rad (DEG (2)), deg2rad (DEG (3)));
        RXYZ = [tree.X tree.Y tree.Z] * RM;
        tree.X = RXYZ (:, 1);
        tree.Y = RXYZ (:, 2);
        tree.Z = RXYZ (:, 3);
    end
end

if regexpi (options, '-s') % show option
    clf; shine; hold on; plot_tree (intree); plot_tree (tree, [1 0 0]);
    HP (1) = plot (1, 1, 'k-'); HP (2) = plot (1, 1, 'r-');
    legend (HP, {'before', 'after'}); set (HP, 'visible', 'off');
    title  ('rotate tree');
    xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
    view (3); grid on; axis image;
end

if (nargout > 0||(isstruct(intree)))
    varargout {1}  = tree; % if output is defined then it becomes the tree
    if nargout > 1 && regexpi (options, '-m3d')
       varargout {2} = angout;
    end
    if nargout > 2 && regexpi (options, '-m3d')
        varargout {3} = cat(2,d,raxis);
    end
    if nargout > 3 && regexpi (options, '-m3d')
        varargout {4} = qual;
    end
else
    trees {intree} = tree; % otherwise add to end of trees cell array
end