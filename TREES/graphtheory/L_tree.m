function L = L_tree(tree)

    L = (sqrt((tree.X-tree.X(idpar_tree(tree))).^2+(tree.Y-tree.Y(idpar_tree(tree))).^2+(tree.Z-tree.Z(idpar_tree(tree))).^2));
