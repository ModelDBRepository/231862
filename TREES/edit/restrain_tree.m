function tree = restrain_tree(tree,maxpl,options)
% restrains ("prunes") tree so that maxpl is not exceeded

% OPTIONS 
% '-i' interpolate termination point to maxpl (default). without this option, the termination point is deleted without substitute

if nargin < 2
    maxpl = 400;
end
if nargin < 3
    options = '-i';
end
PL = Pvec_tree(tree); %get path lengths

if any(PL > maxpl)
    if ~isempty(strfind(options,'-i'))
        idpar = idpar_tree(tree); % get parent indices
        ind = PL > maxpl & PL(idpar) > maxpl; % delete all nodes whose parent nodes are too far away from soma, too
        
        tree = delete_tree(tree,ind,'-r');
        
        % for the rest, make them being as far away as possible (maxpl)
        % without changing direction
        idpar = idpar_tree(tree); % get parent indices of cutted tree
        PL = Pvec_tree(tree); %get path lengths of cutted tree
        ind = PL > maxpl;
        direction = dir_tree(tree,'-n');            % get normed direction vectors
        
        tree.X(ind) = tree.X(idpar(ind)) + direction(ind,1) .* (maxpl-PL(idpar(ind))); %substract path length from parent node and multiply by direction to have point farthest away
        tree.Y(ind) = tree.Y(idpar(ind)) + direction(ind,2) .* (maxpl-PL(idpar(ind)));
        tree.Z(ind) = tree.Z(idpar(ind)) + direction(ind,3) .* (maxpl-PL(idpar(ind)));
    else
        tree = delete_tree(tree,PL > maxpl,'-r'); % delete all nodes which are farther away as maxpl
        
    end
end