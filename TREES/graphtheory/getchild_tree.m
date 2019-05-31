function child = getchild_tree(tree,inodes,options)
if nargin < 3
    options = '';
end
if nargin < 2 || isempty(inodes)
    inodes = 1:numel(tree.X);
end
child = NaN(numel(inodes),2);
[row,col] = find(tree.dA(:,inodes));

for n = 1:numel(inodes)
    child(n,1:sum(col==n)) = row(col==n)';
end

if ~isempty(strfind(options,'-1'))
   child = child(:,1); 
end