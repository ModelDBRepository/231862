function [tree, ind] = bleb_tree(tree, dist)
% adds a bleb (thicker diameter) at a certain distance into the axon
if nargin < 2 || isempty(dist)
    dist = 30;
end

rind = find(~cellfun(@isempty,strfind(tree.rnames,'axon')));

len = Pvec_tree(tree);


bleb = find(len > dist & any(repmat(tree.R,1,numel(rind)) == repmat(rind,numel(tree.R),1),2),1,'first');
ind = getchild_tree(tree,bleb,'-1');
tree = delete_tree(tree,sub_tree(tree,ind),'-r');

dir = dir_tree(tree,'-n');

[tree, ind] = insert_tree(tree,[1 , tree.R(bleb) , tree.X(bleb)+2*dir(bleb,1) , tree.Y(bleb)+2*dir(bleb,2) , tree.Z(bleb)+2*dir(bleb,3) , 2 , bleb]); %insert the bleb (2µm diam)
tree.name = strcat(tree.name,'_bleb');

