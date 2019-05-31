function tree = interpd_tree(tree,ind)

PL = Pvec_tree(tree);
ipar = ipar_tree(tree);

if any(ipar(ind(1),:)== ind(2))
    ind = ipar(ind(1),1:find(ipar(ind(1),:)== ind(2)));
elseif any(ipar(ind(2),:)== ind(1))
    ind = ipar(ind(2),1:find(ipar(ind(2),:)== ind(1)));
else
    errordlg('Indices do not lie on the same path to the root')
    return
end

m = (tree.D(ind(end))-tree.D(ind(1)))/(PL(ind(end))-PL(ind(1)));

tree.D(ind) = m*(PL(ind)-PL(ind(1))) + tree.D(ind(1));
