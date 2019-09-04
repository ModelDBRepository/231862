% script to make region names of trees stereotypic and to add a soma and an axon


[tree,treename,treepath]=load_tree;
if iscell(tree{1})
    tree = tree{1};
end


for t = 1:numel(tree)
    if ~any(strcmp(tree{t}.rnames,'adendIML')) && any(strcmp(tree{t}.rnames,'IML'))
        tree{t}.rnames(strcmp(tree{t}.rnames,'IML')) = {'adendIML'};
    end
    if ~any(strcmp(tree{t}.rnames,'adendMML')) && any(strcmp(tree{t}.rnames,'MML'))
        tree{t}.rnames(strcmp(tree{t}.rnames,'MML')) = {'adendMML'};
    end
    if ~any(strcmp(tree{t}.rnames,'adendOML')) && any(strcmp(tree{t}.rnames,'OML'))
        tree{t}.rnames(strcmp(tree{t}.rnames,'OML')) = {'adendOML'};
    end
    if any(strcmp(tree{t}.rnames,'outside')) || any(strcmp(tree{t}.rnames,'OMLoutside'))
        if ~any(strcmp(tree{t}.rnames,'adendOML'))
            tree{t}.rnames(strcmp(tree{t}.rnames,'outside') | strcmp(tree{t}.rnames,'OMLoutside')) = {'adendOML'};
        end
        ind = find(strcmp(tree{t}.rnames,'outside') | strcmp(tree{t}.rnames,'OMLoutside'));
        tree{t}.R(tree{t}.R == ind) = find(strcmp(tree{t}.rnames,'adendOML'));
        for in = ind+1:numel(tree{t}.rnames)
            tree{t}.R(tree{t}.R == in) = in-1;
        end
        tree{t}.rnames(ind) = [];
    end
    if ~any(strcmp(tree{t}.rnames,'soma'))
        ind = tree{t}.R == find(strcmp(tree{t}.rnames,'SGCL')) | tree{t}.R == find(strcmp(tree{t}.rnames,'GCL'));
        [~,ind2] = max(diff(tree{t}.D(ind),2));
        ind(ind2+1:end) = false;
        tree{t}.R(ind) = numel(tree{t}.rnames)+1;
        tree{t}.rnames = [tree{t}.rnames,'soma'];
%         PL = Pvec_tree(tree{t});
%         type = fittype(sprintf('%g * cos (pi * x / (l / 2)) + %g',tree{t}.D(1)/2,tree{t}.D(1)/2));
%         mycurve=fit(PL(ind(1:ind2)),tree{t}.D(ind(1:ind2)),type,'StartPoint',30);
    elseif sum(tree{t}.R == find(strcmp(tree{t}.rnames,'soma'))) == 1 % soma only defined as one point..
        tree{t} = soma_tree(tree{t},21,21,'-r');
    end
    tree{t} = rot_tree(tree{t},[],'-m3DY-al');
    idpar = idpar_tree(tree{t});
    ind = find(tree{t}.R == find(strcmp(tree{t}.rnames,'adendIML')) & tree{t}.R(idpar) == find(strcmp(tree{t}.rnames,'GCL')));
    if isempty(ind)
        ind = find(tree{t}.R == find(strcmp(tree{t}.rnames,'adendIML')) & tree{t}.R(idpar) == find(strcmp(tree{t}.rnames,'soma')));
    end
    if ~isempty(ind)
        tree{t} = tran_tree(tree{t},[-mean(tree{t}.X(ind)) 0 0]);
    else
       error('Error')
    end
    tree{t} = cap_tree(tree{t},'-i-a');
    
end



save_tree(tree,fullfile(treepath,[treename(1:end-4),'_axon.mtr']))