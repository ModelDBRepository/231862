function tree = calculate_loads(neuron,tree)
% similar as in Hay et al 2013

neuron.params.nseg = 3;  % need to change from lambda calculation to standard nseg definition because of high Ra

ostruct.show = 0;
ostruct.capacitance = 0;

% neuron = [];
% for t = 1:numel(tree)
%     neuron.mech{t}.all.pas = struct('g',0.0001,'e',-80,'Ra',200);
% end
%% get total Rin measured at soma
% neuron = manipulate_Ra(neuron,1,'axon');  %macht kaum unterschied...
g_totalsoma = 1./t2n_passTests(neuron,tree,'',ostruct); % unit in µS  ( 1/MOhm)

%% get only soma Rin
neuron = manipulate_Ra(neuron,1,'all');
neuron = manipulate_Ra(neuron,0,'soma');
g_soma = 1./t2n_passTests(neuron,tree,'',ostruct); % unit in µS  ( 1/MOhm)
neuron = manipulate_Ra(neuron,1,'soma');
neuron = manipulate_Ra(neuron,0,'all');

%% get total Rin measured at AIS
for t = 1: numel(tree)
    PL = Pvec_tree(tree{t});
    ind = find(tree{t}.R == find(strcmp(tree{t}.rnames,'axonh')));
    [~,ind2] = min(abs(PL(ind)-mean(PL(ind))));
    ostruct.recordnode(t) = ind(ind2);%cellfun(@(x) find(x.R == find(strcmp(x.rnames,'axonh')),1,'last'),tree); % get last node of AIS
end
ostruct.stimnode = ostruct.recordnode;
g_totalAIS = 1./t2n_passTests(neuron,tree,'',ostruct); % unit in µS  ( 1/MOhm)

%% get only AIS Rin
neuron = manipulate_Ra(neuron,1,'all');
neuron = manipulate_Ra(neuron,0,'axonh');
g_AIS = 1./t2n_passTests(neuron,tree,'',ostruct); % unit in µS  ( 1/MOhm)
%%
Rho_soma = g_totalsoma./g_soma -1;
Rho_AIS = g_totalAIS./g_AIS -1;

for t = 1:numel(tree)
    tree{t}.Rho_soma = Rho_soma(t);
    tree{t}.Rho_AIS = Rho_AIS(t);
end