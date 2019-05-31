function [tree,idpar] = cap_tree(tree,options)
% options:
% -i1 : 1 node each 1 µm (DEFAULT)
% -a : also add axon
% INACT-so : define soma

if nargin < 2
    options = '-i1';
end
dir = dir_tree(tree,'-n');

   
width = tree.D(1);
fun1 = @(x,y) real(sqrt(x.^2-2*y.^2));
if ~isempty(strfind(options,'-a'))
    length = 1350;
    axondiam = normrnd(0.45,0);%.045);
    scale = normrnd(-0.2,0.02);
    fun2 = @(x,y) 0.8*(x-axondiam)*exp(scale * y)+axondiam;
else
    length = width;
    fun2 = fun1;
end

% if ~isempty(strfind(options,'-so'))  % give root soma region
%     if ~any(strcmp(tree.rnames,'soma'))
%        tree.rnames = [tree.rnames,'soma'];
%        tree.R(1) = numel(tree.rnames);
%    end
% end

if ~isempty(strfind(options,'-i'))  % std is 1 node each µm
    ind = cell2mat(textscan(options,'-i%f'));
    if isempty(ind)
        ind = 1;
    end
    lin = 0:ind:floor(length);
else
    lin = linspace(0,width,ceil(width/mean(len_tree(tree))));
    if ~isempty(strfind(options,'-a'))
        lin(end+1) = 1350;
    end
    
end

idpar = 1;
for l = lin(2:end)
%     d = sqrt(width.^2-lin(l).^2);
    if  fun1(width,l) < fun2(width,l)%tree.D(idpar(end)) <= width/2
        d = fun2(width,l);
    else
        d = fun1(width,l);
        il = l;
    end
    if d > 0 
        [tree,idpar(end+1)] = insert_tree(tree,[0,1,[tree.X(1) tree.Y(1) tree.Z(1)] - l*dir(2,:),d,idpar(end)],'nix');
    else
       warning('diameter was zero or less. skipped') 
    end
end
tree.R(idpar) = tree.R(idpar(1)); % give region same as root

pl = Pvec_tree(tree);

% if ~isempty(strfind(options,'-so'))  % give root soma region
%     tree.R(pl <= width) = tree.R(idpar(1)); % also give nodes in other directions the same region (with distance equal to size of cap that was created)
% end
if ~isempty(strfind(options,'-a'))
   if ~any(strcmp(tree.rnames,'axon'))
       tree.rnames = [tree.rnames,'axon'];
   end
   if ~any(strcmp(tree.rnames,'axonh'))
       tree.rnames = [tree.rnames,'axonh'];
   end
   
   indAIS = pl(idpar)>=il & pl(idpar)< 3* width;
   indaxon = pl(idpar)> 3*width;
   tree.R(idpar(indAIS)) = find(strcmp(tree.rnames,'axonh'));
   tree.R(idpar(indaxon)) = find(strcmp(tree.rnames,'axon'));
end

