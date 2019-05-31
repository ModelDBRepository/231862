% PRUNE_TREE   Prunes short branches.
% (trees package)
%
% tree = prune_tree (intree, radius, region, options)
% -------------------------------------------
%
% Cleans tree of terminating segment smaller than radius. If two sub
% branches are smaller, than the shorter one is deleted first
%
% Input
% -----
% - intree::integer:index of tree in trees or structured tree
% - radius::value: scaling factor for radius delimiter  {DEFAULT: no scaling == 1}
% - options::string: {DEFAULT '-w'}
%     '-s' : show
%     '-w' : waitbar
%
% Output
% ------
% if no output is declared the trees are added in trees
% - tree:: structured output tree
%
% Example
% -------
% prune_tree (sample_tree, 10, '-s')
%
% See also quaddiameter_tree RST_tree
% Uses sort_tree idpar_tree ver_tree
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function [pruned_tree,count,delind] = prune_tree (intree, radius, region, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length(trees); % {DEFAULT tree: last tree in trees cell array}
end;
ver_tree (intree); % verify that input is a tree structure

if (nargin<2)||isempty(radius),
    radius = 10;
end
if nargin < 3 || isempty(region)
    region = unique(intree.R); % all regions of tree
end
if size(region,1) > size(region,2)
    region = region';
end
if (nargin<4)||isempty(options),
    options = '-w';
end

% tree  = sort_tree (intree, '-LO'); % sort tree to be BCT conform, heavy parts left
tree=intree;

typeN = (ones (1, size (tree.dA, 1)) * tree.dA)';  % 0 = TP, 1 = CP, 2 = BP
count = 0 ;
delind = [];  %indices of dendrites to delete
% while 1
%     strahler = strahler_tree(tree);
    Segs = dissect_tree(tree,'-r');
    PL = Pvec_tree(tree);
    ipar = ipar_tree(tree);
%     Tsegs = segs(strahler(segs(:,1)) >= 2 & strahler(segs(:,2)) == 1,:);
    Segs(~any(repmat(tree.R(Segs(:,2)),1,numel(region)) == repmat(region,size(Segs,1),1),2),:) = [];  % delete all segments that are not in regions to be pruned
    while 1
        tsegs = find(typeN(Segs(:,2))==0 & typeN(Segs(:,1))==2); % get indices of all terminal segments
        [len, ind] = min(PL(Segs(tsegs,2))-PL(Segs(tsegs,1))); % find smallest terminal segment
        if len <= radius
            thisBP = Segs(tsegs(ind),1);
%             ind2 = find (typeN (1 : Segs(tsegs(ind),2) - 1)==2 , 1, 'last') : Segs(tsegs(ind),2);  % find last BP between soma and Tseg end..from there to the end...
            delind = cat(2,delind,ipar(Segs(tsegs(ind),2),1:find(ipar(Segs(tsegs(ind),2),:)==thisBP)-1));  % remember this terminal segment for deletion (without BP)
            count = count + 1;
            % this part now handles the virtual deletion of the branch
            ind2 = setdiff(find(Segs(:,1)==thisBP),tsegs(ind)); % find the other subbranch
            Segs(Segs(:,2)==thisBP,2) = Segs(ind2,2); % put that subbranch segment together with its parent segment
            Segs([ind2,tsegs(ind)],:) = [];  % delete original subbranch entry and the "deleted branch"
            typeN(thisBP) = typeN(thisBP) -1;  % the BP is now a CP
        else
            break
        end
    end
    if count > 0
        if ~isempty(strfind (options, '-s'))
            pruned_tree = delete_tree(tree,delind,'-r-s');
        else
            pruned_tree = delete_tree(tree,delind,'-r');
        end
    else
        pruned_tree = tree;
    end
    
%     [len, ind] = min(PL(Tsegs(:,2))-PL(Tsegs(:,1)));
%     if len <= radius
%         tree = delete_tree(tree,find (abs (typeN (1 : Tsegs(ind,2) - 1) - 1), 1, 'last') + 1 : Tsegs(ind,2),'-r');
%         typeN = (ones (1, size (tree.dA, 1)) * tree.dA)';
%         count = count + 1;
%     else
%         break
%     end
% end

% if strfind (options, '-s'),
%     clf; hold on;
%     col = zeros(numel(intree.X),3);
%     col(deleted,1) = 1;
%     plot_tree (intree, col);
%     %     plot_tree (intree, [], [], [], [], '-3l');
%     %     plot_tree (tree, [1 0 0], [], [], [], '-3l');
%     HP (1) = plot (1, 1, 'k-'); HP (2) = plot (1, 1, 'r-');
%     legend (HP, {'before', 'after'});
%     set (HP, 'visible', 'off');
% %     title  ('clean tree');
%     xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
% %     view (2); grid on; axis image;
% end

if (nargout == 0) && ~(isstruct(intree))
    trees {intree} = tree; % otherwise the orginal tree in trees is replaced
end

