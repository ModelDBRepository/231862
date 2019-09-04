function [tree,neuron,treeFilename] = GC_initModel(ostruct)

if ~isfield( ostruct,'reducecells')
    ostruct.reducecells = 0;
end
if ~isfield( ostruct,'newborn')
    ostruct.newborn = 0;
end
if ~isfield( ostruct,'forcecalcload')
    ostruct.forcecalcload = 0;
end

neuron.params.celsius = 24;   % temperature
neuron.params.prerun = 400;   % large-dt prerun to let the system equilibrate
neuron.params.v_init = -90;  % initial membrane voltage
neuron.params.dt = 1;       % standard time step (is automatically changed to a smaller dt in most simulations below)
neuron.params.nseg = 'dlambda';  % number of segments, can be constant or 'dlambda' to adjust it according to the d-lambda rule
neuron.params.v_init = -80;

AHflag = false;
if ostruct.vmodel == 0  % only channels active at subthreshold
    mechoptions = '-p';
elseif ostruct.vmodel > 0
    mechoptions = '-a-p-n';  % active model as published in Beining et al 2017
else
    mechoptions = '-o-a-p';  % old AH99 model
    AHflag = true;
end
if ostruct.newborn
    mechoptions = strcat(mechoptions,'-y');
end
if ostruct.ratadjust
    mechoptions = strcat(mechoptions,'-ra');
end

if ostruct.ratadjust
    str = '_ratadjust';
else
    str = '';
end
if ~isnan(ostruct.usemorph)
    
    switch ostruct.usemorph
        case 1    % SH07 cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile(pwd,'morphos','SH_07_all_repairedandsomaAIS_MLyzed_loadadjusted.mtr'),'file')
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','SH_07_all_repairedandsomaAIS_MLyzed_loadadjusted.mtr'));
            else
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','SH_07_all_repairedandsomaAIS_MLyzed.mtr'));
            end
            
            tname = 'SH07all2';
        case 2   % synth adult mouse cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile(pwd,'morphos','mouse_AAVart_old_pruned_axon_loadadjusted.mtr'),'file')
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','mouse_AAVart_old_pruned_axon_loadadjusted.mtr'));
            else
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','mouse_AAVart_old_pruned_axon.mtr'));
            end
            tname = 'mouse_matGC_art';
            neuron.params.exchfolder = 't2nexchange_aGCmorphsim';
            
        case 3  % synth young mouse cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile(pwd,'morphos','mouse_RVart_pruned_axon_loadadjusted.mtr'),'file')
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','mouse_RVart_pruned_axon_loadadjusted.mtr'));
            else
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','mouse_RVart_pruned_axon.mtr'));
            end
            tname = 'mouse_abGC_art';
            neuron.params.exchfolder = 't2nexchange_aGCmorphsim2';
        case 4  %  adult rat cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile(pwd,'morphos',sprintf('Beining_AAV_contra_MLyzed_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos',sprintf('Beining_AAV_contra_MLyzed_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','Beining_AAV_contra_MLyzed_axon.mtr'));
            end
            tname = 'rat_mGC_Beining';
            neuron.params.exchfolder = 't2nexchange_aGCmorphsim6';
        case 5 % synth adult rat cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile(pwd,'morphos',sprintf('rat_AAVart_old_pruned_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos',sprintf('rat_AAVart_old_pruned_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','rat_AAVart_old_pruned_axon.mtr'));
            end
            tname = 'rat_matGC_art';
            neuron.params.exchfolder = 't2nexchange_aGCmorphsim3';
        case 6  % synth young rat cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile(pwd,'morphos',sprintf('rat_RVart_pruned_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos',sprintf('rat_RVart_pruned_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','rat_RVart_pruned_axon.mtr'));
            end
            tname = 'rat_abGC_art';
            neuron.params.exchfolder = 't2nexchange_aGCmorphsim4';
        case 7  % Claiborne rat cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile(pwd,'morphos',sprintf('Claiborne_male_MLyzed_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos',sprintf('Claiborne_male_MLyzed_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treeFilename,treepath]=load_tree(fullfile(pwd,'morphos','Claiborne_male_MLyzed_axon.mtr'));
            end
            tname = 'rat_mGC_Claiborne';
            neuron.params.exchfolder = 't2nexchange_aGCmorphsim5';
            
    end
    
    neuron.experiment = tname;
    if AHflag
        neuron.experiment = strcat(neuron.experiment,'_AH99');
    end
        ntree = min(numel(tree),15);
        tree=tree(1:ntree);
        if ostruct.usecol
            if any(ostruct.usemorph == [1,4]) % reconstructed morphologies
                colors = colorme(ntree,'-grb');
            else
                colors = colorme(ntree,'-grg');
            end
        else
            colors = colorme(ntree);
        end
        for t = 1:numel(tree)
            tree{t}.col = colors(t);
        end
else
    [tree,treeFilename,treepath]=load_tree;
    if isempty(tree)
        return
    end

    colors = colorme(numel(tree));
    for t = 1:numel(tree)
        tree{t}.col = colors(t);
    end
    tname = treeFilename(1:end-4);
end

if ~all(cellfun(@(x) isfield(x,'NID'),tree)) || ~all(cellfun(@(x) exist(fullfile(pwd,'morphos','hocs',[x.NID,'.hoc']),'file'),tree))
    answer = questdlg('Caution! Not all of your trees have been transformed for NEURON yet! Transforming now..','Transform trees','OK','Cancel','OK');
    if strcmp(answer,'OK')
        tree = t2n_writeTrees(tree,tname,fullfile(treepath,treeFilename));
    end
end
if isempty(tree)
    return
end
if isstruct(tree)
    tree={tree};
elseif iscell(tree{1})
    tree=tree{1};
end

cd(pwd)

for t = 1:numel(tree)
    tree{t} = sort_tree(tree{t},'-LO');
    
    neuron.mech{t} = [];
    if ~AHflag && ostruct.scalespines  % spine scaling ignored in AH99 as already taken into account in the biophys
        neuron.mech{t} =  t2n_catStruct(GC_biophys(mechoptions),GC_spinedensity(ostruct.scalespines*0.9));
    else
        neuron.mech{t} =  t2n_catStruct(GC_biophys(mechoptions));
    end
    if ~AHflag || ostruct.changeAHion
        neuron = t2n_setionconcentration(neuron,'Mongiat');
    end
    if ~isfield(tree{t},'col')
        tree{t}.col{1} = rand(1,3);
    end
    if ostruct.noise ~= 0
        neuron.pp{t}.InGauss = struct('node',1,'mean',0.01,'stdev',0.01,'del',0,'dur',1e9);
    end
end

if ostruct.scalespines
    neuron = t2n_implicitSpines(neuron);
end
%%% This is the Hay 2013 implementation of adjusting soma and AIS
%%% conductance according to dendritic morphology
if ostruct.adjustloads
    if any(cellfun(@(x) ~isfield(x,'Rho_soma') | ~isfield(x,'Rho_AIS'),tree)) || ostruct.forcecalcload
        tree = calculate_loads(neuron,tree);
        save_tree(tree,fullfile(treepath,[treeFilename(1:end-4),sprintf('_loadadjusted%s.mtr',str)]));
    end
    if ostruct.usemorph >= 4      %rat
        neuron = adjust_loads(neuron,tree,'r',ostruct);
    else        %mouse
        neuron = adjust_loads(neuron,tree,'m');
    end
end

if ostruct.reducecells
    if ostruct.usemorph ~= 1
        tree=tree((1:3)+2);
        neuron.mech = neuron.mech((1:3)+2);
    else
        tree=tree(3);
        neuron.mech = neuron.mech(3);
    end
    neuron.experiment = strcat(neuron.experiment,'_reduceNs');
end
if ostruct.ratadjust
    neuron.experiment = strcat(neuron.experiment,'_ratadjust');
end

treeFilename = fullfile(treepath,treeFilename);
fprintf('Model initialized...%s\n',neuron.experiment)

