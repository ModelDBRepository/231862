%% comments

% This is the script containing all simulation protocols for the figures of
% the paper. It can either be run completely, or single sections ( =
% simulation protocols) can be run separately, when the first section
% (initialization) has been run before

%% initialize trees and mechanisms

%********* folders
targetfolder_results = 'C:/GCModel/results';  % the folder where the graphs and pictures are saved
targetfolder_data = 'C:/GCModel/simdata';  % the folder where simulated data is saved (to avoid resimulating if not necessary)
% ***********

ostruct = struct('plot','auto','show',3,'legend',0,'marker','o','sem',1,'FontSize',10,'FontType','Arial','barwidth',0.3,'figurewidth',8,'figureheight',5,'LineWidth',1,'grid',0,'priority','plot','ticklength',0.015);  % some options that are used when something is plotted
ostruct.usecol = 1;  % 0 = pseudorandom colors for each simulated cell in graph, 1 = green or blue grading


% change biophysical model here
ostruct.vmodel = 1; % 0 = passive model, > 0 = active model, everything else (e.g. NaN) = old AH99 model
ostruct.changeAHion = 0;  % only important when using the AH99 model. Boolean to decide if standard AH99 ion reversal potentials are used (0) or if they are adjusted to the experiments (1)

% change morphologies here
ostruct.usemorph = 1;  % 1 = all SH07, 2= synth mouseMat, 3= synth mouseYoung 4= Beining (2016) AAV rat, 5 = synth ratOld 6= synth ratYoung 7 = Claiborne,
ostruct.newborn = 0;  % 0 = adult GC model, 1 = young abGC model

% more parameters
ostruct.reducecells = 0;  % reduce number of cells for faster simulation (e.g. for testing)
ostruct.scalespines = 1;  % scaling of g_pas and cm to implicitly model spines. is ignored in AH99 model, as this is already taken into account in the biophys model!
ostruct.adjustloads = 0;  % the Hay et al 2013 implementation of adjust dendritic loads to reduce variability between cells (not used in publication)
ostruct.noise = 0;       % add noise to the membrane voltage by injecting gaussian noise current to the soma (not working with variable dt / cvode)

% do not change from here
if ostruct.usemorph >= 4
    ostruct.ratadjust = 1;  % adjust Kir channel in rats
else
    ostruct.ratadjust = 0;
end

%*****************************
if ostruct.newborn
    ostruct.scalespines = 0.3 * ostruct.scalespines;  % means g_pas and cm are only scaled by 30% of the original spine densities due to reduced spine density in young abGCs
end
if ~exist(fullfile(pwd,'GC_initModel.m'),'file')
    error('It seems your current folder is not the GC model folder. Please change it')
end
if ~exist(targetfolder_data,'file')
    mkdir(targetfolder_data)
end
if ~exist(targetfolder_results,'file')
    mkdir(targetfolder_results)
end
[tree,neuron,treename] = GC_initModel(ostruct);  % initialize the model by loading the morphologies and setting the biophysical parameters
if ostruct.newborn
    %     % this is the good FI,deep fAHP one...
    %     ostruct.channelblock = {'Kir21','Kv42','na8st','BK','Cav13','Kv21'};%,'Kv21','Kv42','Kv14','Kv34','na8st','Cav13'};%'Cav22'};     %{'Kv42','Kir21','Kv14','Kv34','pas','Kv21','na8st','Kv723'};%,'na8st','Kv723','na8st','Kv21'};%'except','Kir21'};%{'Kir','Kv42','HCN'};%'Kir','Kv42','Kv14','HCN'}; % Kir , SK etc
    %     ostruct.blockamount = [73,50,25,60,50,75];%68 kir %,100,80,80,50,55,100];%  Kv42 50 Kv21 70
    %      ostruct.specify = {'','','','','',''};
    
    %      % this is the used version
    oldexpname = neuron.experiment;
    neuron = t2n_blockchannel(neuron,{'Kir21','Kv42','na8st','BK','Cav13','Kv21','Kv723','BK'},[73,50,25,40,50,50,50,100],[],{'','','','gakbar','','','','gabkbar'});
    neuron.experiment = strcat(oldexpname,'_newborn');
    disp('young abGC version')
end
neuron_orig = neuron;

%% plot trees, Figure 1
% ostruct.show = 2;
% ostruct.savename = 'Fig1_Trees';
% t2n_plotTrees(tree,targetfolder_results,[],ostruct) % '-s'
if ostruct.usemorph < 4  % mouse experiments
    %% Mongiat IV + Ba, simulate the I-V curve with and without blocking Kir channels with Barium, Figure 2 & 6
    neuron = neuron_orig;
    
    holding_voltage = -80 - 12.1; % mV, LJP corrected
    neuron.params.cvode = 1;  % boolean if dt is constant (0) or variable (1)
    neuron.params.dt = 0.25;  % this is ignored if cvode = 1
    vstepsModel =  (-130:5:-40) - 12.1; % LJP corrected
    dur = [105 100 105];
    
    ostruct.subtract_hv = 1; % boolean subtract holding voltage current
    ostruct.single = 0;      % 1 = show I-V of each cell
    ostruct.figureheight = 4;
    ostruct.figurewidth = 6;
    
    %
    t2n_voltSteps(neuron,tree,vstepsModel,dur,holding_voltage,targetfolder_data);
    %
    if ~ostruct.newborn
        ostruct.dataset =3;  % 1 = old mature dataset, 2 = old young dataset, 3 = new mature dataset, 4 = new young dataset, 5 = new mature BaCl dataset, 6 = new young BaCl dataset
        ostruct.savename = sprintf('Fig2-IV_dyn-%s',neuron.experiment);
        ostruct.handles = t2n_plotVoltSteps(targetfolder_data,neuron,[],ostruct.subtract_hv);
        aGC_plotVoltStepsExp(t2n_catName(targetfolder_data,'Exp_VoltSteps',neuron.experiment,'.mat'),targetfolder_results,ostruct);
        ostruct.handles = [];
        savename = sprintf('Fig2-IV+Ba-%s',neuron.experiment);
    else
        ostruct.dataset = 2.28;
        savename = sprintf('Fig6-IV+Ba_young-%s',neuron.experiment);
    end
    
    ostruct.handles = t2n_IVplot(targetfolder_data,neuron_orig,ostruct);
    aGC_IVplotExp(ostruct,t2n_catName(targetfolder_data,'Exp_VoltSteps',neuron_orig.experiment,'.mat'))
    
    neuron = t2n_blockchannel(neuron,{'Kir21','pas'},[99 30]);
    t2n_voltSteps(neuron,tree,vstepsModel,dur,holding_voltage,targetfolder_data);
    
    t2n_IVplot(targetfolder_data,neuron,ostruct)
    if ~ostruct.newborn % plot mongiat data
        ostruct.dataset =5;
        aGC_IVplotExp(ostruct,t2n_catName(targetfolder_data,'Exp_VoltSteps',neuron_orig.experiment,'.mat'))
    end
    ostruct.dataset =0; % plot literature data (arrows)
    aGC_IVplotExp(ostruct,t2n_catName(targetfolder_data,'Exp_VoltSteps',neuron_orig.experiment,'.mat'))
    
    ostruct.handles = [];
    xlim([-140 -60])
    if ostruct.newborn
        ylim([-150 50])
    else
        ylim([-400 100])
    end
    FontResizer
    tprint(fullfile(targetfolder_results,savename),'-pdf');
    
    %% Mongiat FI + Ba, simulate the F-I relationship with and without blocking Kir by application of Barium, Figure 3 & 7
    neuron = neuron_orig;
    neuron.experiment = strcat(neuron.experiment,'_Fig37');
    ostruct.holding_voltage = -80; % mV % -80
    % exp_iclamp = load_ephys(ostruct.dataset,'CClamp');  % alternatively get holding voltage from experiment
    % ostruct.holding_voltage = mean(mean(exp_iclamp(:,:,1)));
    
    ostruct.figureheight = 4;
    ostruct.figurewidth = 6;
    ostruct.amp = (0:5:120)/1000;% current steps in nA which are to be simulated
    neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
    ostruct.coarse = 0.5;  % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
    
    if ostruct.newborn
        ostruct.dataset = 2.28;
        ostruct.savename = sprintf('Fig7_Final-%s',neuron.experiment);
        ostruct.savename2 = sprintf('Fig7-FI-%s_young',neuron.experiment);
        ostruct.savename3 = sprintf('Fig7-FI+Ba-%s_young',neuron.experiment);
        steps = [10, 50]/1000;  % current step [pA] of which voltage trace is plotted
    else
        ostruct.dataset = 3;
        ostruct.savename = sprintf('Fig3_Final-%s',neuron.experiment);
        ostruct.savename2 = sprintf('Fig3-FI-%s',neuron.experiment);
        ostruct.savename3 = sprintf('Fig3-FI+Ba-%s',neuron.experiment);
        if isnan(ostruct.vmodel) % AH99 model. This is different to our model because no spiking occurred at 75 pA in the AH99 GC model
            steps = [30, 90]/1000;  % current step [pA] of which voltage trace is plotted
        else
            steps = [30, 75]/1000;  % current step [pA] of which voltage trace is plotted
        end
    end
    
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    
    ostruct.handles = figure;
    t2n_plotCurrSteps(targetfolder_data,neuron,steps);
    if ~ostruct.reducecells && numel(ostruct.amp) == 25
        vec = true(numel(tree)*2,1);
        if ostruct.newborn
            vec(cat(2,[1,3,7],[1,3,7]+1)) = false;
%             delete(ostruct.handles.Children.Children(reshape(logical(repmat([1,0,1,1,0,0,1,1],2,1)),16,1))) % only keep 3 of the 8 morphs for visibility
        else
            switch ostruct.usemorph
                case 1
                    vec(cat(2,[9,11,15],[9,11,15]+1)) = false;
                case 2
                    vec(cat(2,[1,27,29],[1,27,29]+1)) = false;
            end
%             delete(ostruct.handles.Children.Children(reshape(logical(repmat([1,1,1,1,0,0,1,0],2,1)),16,1))) % only keep 3 of the 8 morphs for visibility
        end
        delete(ostruct.handles.Children.Children(vec))
    end
    FontResizer
    FigureResizer(ostruct.figureheight,ostruct.figurewidth)
    xlim([0 350])
    ylim([-85 -30])
    if ostruct.newborn
        tprint(fullfile(targetfolder_results,strcat('Fig.7-SpikingModel_',neuron.experiment)),'-pdf');
    else
        tprint(fullfile(targetfolder_results,strcat('Fig.3-SpikingModel_',neuron.experiment)),'-pdf');
    end
    ostruct.handles = figure;
    t2n_plotCurrSteps(targetfolder_data,neuron);
    aGC_plotCurrStepsExp(targetfolder_data,targetfolder_results,neuron,ostruct,steps)
    ostruct.handles = [];
    % % this part is only used if experimental data and convex hulls of data
    % % should be plotted (for publication)
    % ostruct.savename = strcat(ostruct.savename,'_onlyexp');
    % [prop1,fig] =     t2n_APprop(targetfolder_data,neuron,ostruct);
    %     ostruct.blur = 1;
    %     ostruct.patch = 1;
    %     figure(fig(6))
    %     [~,fig2] = plot2dens(gca,[2 2 1],ostruct);
    %     figure(fig2(4))
    %     tprint(fullfile(targetfolder_results,strcat(ostruct.savename,'-','APvthreshVSAPwidth_densitypatch')),'-pdf');
    %     figure(fig(7))
    %     [~,fig2] = plot2dens(gca,[20 0.2 1],ostruct);
    %     figure(fig2(4))
    %     tprint(fullfile(targetfolder_results,strcat(ostruct.savename,'-','APampVSfahp_densitypatch')),'-pdf');
    % ostruct.savename = ostruct.savename(1:end-8);
    
    [prop2,fig,figname] =     t2n_APprop(targetfolder_data,neuron,0.09,tree,ostruct.dataset);
    xl = repmat([0.5 5.5],5,1);
    if ostruct.newborn
        yl = [0,-60,-80,0,0;4,-20,-40,100,150]';
    else
        if isnan(ostruct.vmodel) %AH99
            yl = [0,-60,-80,0,0;6,-30,-40,100,150]';
        else
            yl = [0,-60,-80,0,0;2,-30,-40,100,150]';
        end
    end
    for f = 1:14
        figure(fig(f))
        if f <= 5
            ylim(yl(f,:))
            xlim(xl(f,:))
            set(gca,'XTick',1:5)
        else
            switch f
                case 6
                    xlim(yl(2,:))
                    if ostruct.newborn
                        ylim([-10 30])
                    else
                        ylim([0 25])
                    end
                case 7
                    xlim(yl(5,:))
                    ylim(yl(1,:))
                case 8
                    xlim([20 70])
                    if ostruct.newborn
                        ylim([-70 -35])
                    else
                        ylim([-70 -45])
                    end
                case 9
                    xlim(xl(1,:))
                    ylim([0 4])
                    set(gca,'XTick',1:5)
            end
        end
        FontResizer
        FigureResizer(ostruct.figureheight,ostruct.figurewidth)
        if isfield(ostruct,'savename')
            tprint(fullfile(targetfolder_results,strcat(ostruct.savename,'-',figname{f})),'-pdf');
        else
            tprint(fullfile(targetfolder_results,strcat(figname{f}),neuron.experiment),'-pdf');
        end
    end
    
    ostruct.handles = [];
    ostruct.savename = ostruct.savename2;
    if ~ostruct.newborn
        ostruct.figurewidth = 4;
    end
    ostruct.handles = t2n_FIplot(targetfolder_data,neuron,ostruct,targetfolder_results);
    aGC_FIplotExp(targetfolder_data,targetfolder_results,neuron,ostruct)

    if ~ostruct.newborn
        neuron = t2n_blockchannel(neuron,{'Kir21','pas'},[99 30]);
        if ostruct.newborn
            ostruct.dataset =0;
        else
            ostruct.dataset =5;
        end
        
        t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
        ostruct.handles = [];
        ostruct.savename = ostruct.savename3;
        ostruct.handles = t2n_FIplot(targetfolder_data,neuron,ostruct,targetfolder_results);
        aGC_FIplotExp(targetfolder_data,targetfolder_results,neuron,ostruct)
        
        ostruct.handles = [];
        ostruct.savename = [];
    end
    
    
    %% Mongiat dV / dt curve,Figure 3 & 7, phase plots, extra simulation needed because a small dt is necessary to correctly assess maximal dV/dt's
    neuron = neuron_orig;
    
    neuron.experiment = strcat(neuron.experiment,'_dV');
    ostruct.amp = [40,65,90,115]/1000; % steps that are simulated plotted [nA]
    ostruct.holding_voltage = -80;
    ostruct.duration = 200;
    neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
    ostruct.coarse = 0;   % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
    ostruct.figurewidth = 6;
    ostruct.figureheight = 4;
    ostruct.ampprop = 0.090;  % plot this current step [nA] in the single figure
    if ostruct.newborn
        ostruct.dataset =2.28;
        ostruct.savename = sprintf('Fig7_-%s',neuron.experiment);
    else
        ostruct.dataset =3;
        ostruct.savename = sprintf('Fig3_-%s',neuron.experiment);
    end
    if ostruct.usemorph >= 4  % = rat morphology
        ostruct.dataset = 0;
        ostruct.savename = sprintf('SupplFigX_-%s',neuron.experiment);
    end
    
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    
    [maxdv,ostruct.handles] = t2n_plotdV(targetfolder_data,neuron,ostruct,targetfolder_results);  % plot the dV plot
    aGC_plotdVExp(targetfolder_results,ostruct);
    
    %% Kv1.1 overexpression, Figure 6, reproduce Kirchheim et al 2013, Kv1.1 overexpression after status epilepticus reduced FI and increases spiking delay
    steps = 70/1000;
    neuron = neuron_orig;
    ostruct.handles = [];
    ostruct.holding_voltage = -80; %mV
    neuron.experiment = strcat(neuron.experiment,'_Fig6');
    ostruct.figureheight = 4;
    ostruct.figurewidth = 6;
    ostruct.duration = 200;
    ostruct.amp = (0:5:120)/1000; % nA
    neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
    ostruct.coarse = 0.5;   % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
    
    ostruct.dataset =3;
    ostruct.savename1 = sprintf('Fig6_Kv11-%s',neuron.experiment);
    ostruct.savename2 = sprintf('Fig6-FI-Kv11-%s',neuron.experiment);
    
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    handles1 = figure;
    t2n_plotCurrSteps(targetfolder_data,neuron,steps);
    ostruct.color = [0 0 1];
    ostruct.handles2 = t2n_FIplot(targetfolder_data,neuron,ostruct,targetfolder_results);
    newcol = colorme(numel(neuron.mech),'-grr');
    tmptree = tree;
    %increase Kv1.1 10fold
    for t = 1:numel(neuron.mech)
        fields = fieldnames(neuron.mech{t});
        for f1 = 1:numel(fields)
            if isfield(neuron.mech{t}.(fields{f1}),'Kv11')
                neuron.mech{t}.(fields{f1}).Kv11.gkbar = neuron.mech{t}.(fields{f1}).Kv11.gkbar * 10;
            end
        end
        tmptree{t}.col = newcol(t);  % give trees different color for the plot
    end
    neuron.experiment = strcat(neuron.experiment,'_Kv11oe');
    t2n_currSteps(neuron,tmptree,targetfolder_data,ostruct)  % do the simulation
    
    ostruct.savename = ostruct.savename1;
    figure(handles1);
    t2n_plotCurrSteps(targetfolder_data,neuron,steps);
    
    ostruct.savename = ostruct.savename2;
    ostruct.handles = ostruct.handles2;
    ostruct.color = [1 0 0];
    figure(ostruct.handles2(1));ylim([0 8]),set(gca,'YTick',0:2:8)
    t2n_FIplot(targetfolder_data,neuron,ostruct,targetfolder_results);
    ostruct.handles = [];
    ostruct.savename = [];
    
    %% Synaptic stimulation, Figure 8
    ostruct.handles = [];
    neuron = neuron_orig;
    type = 'temporal'; % test, theta-burst stimulation, regular input, temporal shift in input, spatial shift in input
    ostruct.synmode = 1; % 1 = no adjustment of synapse number, 2 = newborns have only half of synapse number, 3 = adjust the number of synapses to get subthreshold response
    % ostruct.figureheight = 3;
    aGC_synstim(neuron,tree,type,ostruct,targetfolder_data)  % do the simulation
    
    ostruct.handles = aGC_synstim_plot(neuron,type,ostruct,targetfolder_data,targetfolder_results,0);
    
else   % to avoid running through rat experiments when just running the script
    %% *********************************  RAT EXPERIMENTS **********************************
    
    %% passive parameter assessment rat
    neuron = neuron_orig;
    ostruct.passtest = 'Std'; % different protocols to measure Rin/tau/capacitance/Vrest (not all can do everything): Mongiat Mongiat2 SH Brenner Std
    [Rin, tau, cap, Vrest] = t2n_passTests(neuron,tree,targetfolder_results,ostruct);
    %% IV rat, generate I-V relationship curve and additinally plot Rin measurements from rat, Figure 4
    neuron = neuron_orig;
    vstepsModel =  -130:5:-40;
    dur = [105 100 105];
    holding_voltage = -80; % mV
    ostruct.subtract_hv = 1; % boolean subtract holding voltage current
    neuron.params.cvode = 1;  % boolean if dt is constant (0) or variable (1)
    ostruct.coarse = 0;   % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
    ostruct.single = 0; % show single data curves instead of mean
    
    t2n_voltSteps(neuron,tree,vstepsModel,dur,holding_voltage,targetfolder_data);  % do the simulation
    
    ostruct.savename = sprintf('Fig4-IV_Mehranfard15-%s',neuron.experiment);
    if ostruct.ratadjust && isempty(strfind(neuron.experiment,'AH99'))
        ostruct.savename = strcat(ostruct.savename,'_ratadjust');
    end
    ostruct.dataset = 0;
    ostruct.handles = t2n_IVplot(targetfolder_data,neuron,ostruct);
    xlim([-125 -60])
    if ostruct.newborn
        ylim([-200 100])
    else
        ylim([-400 200])
    end
    FontResizer
    aGC_IVplotExp(ostruct,t2n_catName(targetfolder_data,'Exp_VoltSteps',neuron.experiment,'.mat'))
    
    tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
    
    %% Rat FI curve, Figure 4, similar to Mehranfard 2015 Plos One
    neuron = neuron_orig;
    if isfield(ostruct,'handles')
        ostruct = rmfield(ostruct,'handles');
    end
    ostruct. handles = [];
    ostruct.amp = (50:50:300)/1000; % nA
    ostruct.duration = 1000;
    neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
    ostruct.coarse = 0.5;   % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
    ostruct.holding_voltage = -80;  % unknown. but without LJP
    
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    
    ostruct.dataset = 0;  % 1 = old mature dataset, 2 = old young dataset, 3 = new mature dataset, 4 = new young dataset, 5 = new mature BaCl dataset, 6 = new young BaCl dataset
    ostruct.savename = sprintf('Fig4_FI_Mehranfahrd15-%s',neuron.experiment);
    if ostruct.ratadjust && isempty(strfind(neuron.experiment,'AH99'))
        ostruct.savename = strcat(ostruct.savename,'_ratadjust');
    end
    ostruct.handles = t2n_FIplot(targetfolder_data,neuron,ostruct,targetfolder_results);
    aGC_FIplotExp(targetfolder_data,targetfolder_results,neuron,ostruct)
    
    ostruct.savename = sprintf('Fig4_Final-%s',neuron.experiment);
    if ostruct.ratadjust && isempty(strfind(neuron.experiment,'AH99'))
        ostruct.savename = strcat(ostruct.savename,'_ratadjust');
    end
    figure; hold all
    t2n_plotCurrSteps(targetfolder_data,neuron,0.1);
    
    %% ISI spike adaptaption Figure 6, similar to MA14 Fig 8
    neuron = neuron_orig;
    %!!!
    neuron.params.celsius = 33;
    %!!!
    ostruct.holding_voltage = -62;
    ostruct.numAP = 6; % 6 spikes per current step
    % ostruct.amp = 200/1000;%nA this value is ignored because of numAP
    neuron.experiment = strcat(sprintf('MA14Fig8_freq%g_%dC_',ostruct.numAP,neuron.params.celsius),neuron.experiment);
    
    if isfield(ostruct,'handles')
        ostruct = rmfield(ostruct,'handles');
    end
    ostruct. handles = [];
    steps = 0.15;%,0.2];
    ostruct.coarse = 0.5;   % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
    ostruct.duration = 100;
    neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
    ostruct.ampprop = 200/1000;
    ostruct.savename = 'Fig6-ISIadapt100ms_MA14';
    if ostruct.ratadjust && isempty(strfind(neuron.experiment,'AH99'))
        ostruct.savename = strcat(ostruct.savename,'_ratadjust');
    end
    ostruct.dataset = 0;  % 1 = old mature dataset, 2 = old young dataset, 3 = new mature dataset, 4 = new young dataset, 5 = new mature BaCl dataset, 6 = new young BaCl dataset
    
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    fig = figure; hold all
    t2n_plotCurrSteps(targetfolder_data,neuron,steps);
    ax = gca; p1 = ax.Children; for n = 1:numel(p1), p1(n).Color = [0 0 0]; end
    props = t2n_APprop(targetfolder_data,neuron,ostruct.ampprop);
    
    ISIadp{1} = NaN(numel(tree),1);
    for t = 1:numel(tree)
        if numel(props.APISI{1,t}) > 1
            ISIadp{1}(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
        end
    end
    
    neuron = t2n_blockchannel(neuron,{'SK2'},100);
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    figure(fig)
    t2n_plotCurrSteps(targetfolder_data,neuron,steps);
    p2 = setdiff(ax.Children,p1); for n = 1:numel(p2), p2(n).Color = [1 0 0]; end
    props = t2n_APprop(targetfolder_data,neuron,ostruct.ampprop);
    
    ISIadp{2} = NaN(numel(tree),1);
    for t = 1:numel(tree)
        if numel(props.APISI{1,t}) > 1
            ISIadp{2}(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
        end
    end
    
    neuron = neuron_orig;
    neuron.experiment = strcat(sprintf('MA14Fig8_freq%g_%dC_',ostruct.numAP,neuron.params.celsius),neuron.experiment);
    neuron = t2n_blockchannel(neuron,{'Kv723'},100);
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    
    figure(fig)
    t2n_plotCurrSteps(targetfolder_data,neuron,steps);
    p3 = setdiff(ax.Children,[p1,p2]); for n = 1:numel(p3), p3(n).Color = [0 0 1]; end
    
    props = t2n_APprop(targetfolder_data,neuron,ostruct.ampprop);
    
    ISIadp{3} = NaN(numel(tree),1);
    for t = 1:numel(tree)
        if numel(props.APISI{1,t}) > 1
            ISIadp{3}(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
        end
    end
    
    dstruct = [];gstruct = [];
    dstruct.data.ISIadapt = [];
    dstruct.strs.ISIadapt.ylab = 'Adaptation index [%]';
    gstruct.group{1} = [ones(1,numel(tree)), repmat(2,1,numel(tree)), repmat(3,1,numel(tree))];
    gstruct.ugroupdef{1} = {'CTRL','SK block','Kv7 block'};
    dstruct.strs.ISIadapt.xstr = gstruct.ugroupdef{1};
    
    for n = 1:3
        dstruct.data.ISIadapt(1,(n-1)*numel(tree)+1:numel(tree)*n) = ISIadp{n}./ISIadp{1}*100;
    end
    
    Masterplotter(dstruct,gstruct,[],ostruct);
    neuron = neuron_orig;
    neuron.experiment = strcat(sprintf('MA14Fig8_freq%g_%dC_',ostruct.numAP,neuron.params.celsius),neuron.experiment);
    tprint(t2n_catName(targetfolder_results,'Fig.6-ISIadapt',neuron.experiment),'-pdf')
    
    neuron.params.celsius = 24;
    ostruct = rmfield(ostruct,'numAP');
    
    %% AP width BK/Kv34 + spike adaptation, Figure 6
    neuron = neuron_orig;
    neuron.experiment = strcat(neuron.experiment,'_APwidth');
    ostruct.amp = [90,250]/1000; % nA  < 10 HZ and höher
    ostruct.duration = 1000;
    neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
    ostruct.coarse = 0.5;   % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
    ostruct.holding_voltage = -80;
    ostruct.ampprop = 90/1000;
    if ostruct.newborn
        ostruct.savename = 'Fig6-APprop_newborn';
    else
        ostruct.savename = 'Fig6-APprop';
    end
    
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    
    ostruct.handles = [];
    prop = cell(2,1);
    for a = 1:2
        prop{a}(1) = t2n_APprop(targetfolder_data,neuron,ostruct.amp(a));
    end
    
    ostruct.handles = figure; hold all
    t2n_plotCurrSteps(targetfolder_data,neuron)
    
    
    for t = numel(tree):-1:1
        if t == 1
            ostruct.handles(1).Children(end-1).Children(t).Color = [0 0 0];
            ostruct.handles(1).Children(end).Children(t).Color = [0 0 0];
        else
            delete(ostruct.handles(1).Children(end-1).Children(t))
            delete(ostruct.handles(1).Children(end).Children(t))
        end
    end
    
    neuron = t2n_blockchannel(neuron,{'Kv34'},100);
    %
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)  % do the simulation
    %
    for a = 1:2
        prop{a}(2) = t2n_APprop(targetfolder_data,neuron,ostruct.amp(a));
    end
    figure(ostruct.handles)
    t2n_plotCurrSteps(targetfolder_data,neuron)
    
    for t = numel(tree):-1:1
        if t == 1
            ostruct.handles(1).Children(end-1).Children(t).Color = [1 0 0];
            ostruct.handles(1).Children(end).Children(t).Color = [1 0 0];
        else
            delete(ostruct.handles(1).Children(end-1).Children(t))
            delete(ostruct.handles(1).Children(end).Children(t))
        end
    end
    
    neuron = neuron_orig;
    neuron.experiment = strcat(neuron.experiment,'_APwidth');
    neuron = t2n_blockchannel(neuron,{'BK'},100);
    %
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)
    %
    for a = 1:2
        prop{a}(3) = t2n_APprop(targetfolder_data,neuron,ostruct.amp(a));
    end
    figure(ostruct.handles)
    t2n_plotCurrSteps(targetfolder_data,neuron)
    
    for t = numel(tree):-1:1
        if t == 1
            ostruct.handles(1).Children(end-1).Children(t).Color = [0 0 1];
            ostruct.handles(1).Children(end).Children(t).Color = [0 0 1];
        else
            delete(ostruct.handles(1).Children(end-1).Children(t))
            delete(ostruct.handles(1).Children(end).Children(t))
        end
    end
    % % SK/Kv7 block
    %
    % %
    dstruct = [];gstruct = [];
    dstruct.data.APwidth = [];
    dstruct.data.APwidth2 = [];
    dstruct.strs.APwidth.xstr = sprintfc('%d pA',ostruct.amp*1000 );
    dstruct.strs.APwidth.ylab = 'Half-amplitude width [ms]';
    dstruct.strs.APwidth2.xstr = sprintfc('%d pA',ostruct.amp*1000 );
    dstruct.strs.APwidth2.ylab = 'Spike width [ms]';
    gstruct.group{1} = [ones(1,numel(tree)), repmat(2,1,numel(tree)), repmat(3,1,numel(tree))];%[4,4,5,5,6,6,ones(1,numel(tree)), repmat(2,1,numel(tree)), repmat(3,1,numel(tree))];
    gstruct.ugroupdef{1} = {'CTRL','Kv3.4 block','BK block'};%,'CTRL data','Kv3.4 block data','BK block data'};
    xth = [1,5];
    for a = 1:2
        for n = 1:3
            dstruct.data.APwidth(a,(n-1)*numel(tree)+1:numel(tree)*n) = cellfun(@(x) x(xth(a)),prop{a}(n).APwidth); % get AP widths of xth AP in all simulations
            dstruct.data.APwidth2(a,(n-1)*numel(tree)+1:numel(tree)*n) = cellfun(@(x) x(xth(a)),prop{a}(n).APwidth2); % get AP widths of xth AP in all simulations
        end
    end
    additionalData = [0.9662-0.0185,0.9662+0.0185,NaN,NaN,1.5748-0.0154,1.5748+0.0154;1.5345+0.1081,1.5345-0.1081,NaN,NaN,2.1089-0.3732,2.1089+0.3732];
    additionalData(:,3:4) = [mean(additionalData(:,1:2),2)*(1.0793-0.0205),mean(additionalData(:,1:2),2)*(1.0793+0.0205)];
    additionalData = cat(2,additionalData(:,1:2),NaN(2,numel(tree)-2),additionalData(:,3:4),NaN(2,numel(tree)-2),additionalData(:,5:6),NaN(2,numel(tree)-2));
    dstruct.data.APwidth = cat(1,additionalData,dstruct.data.APwidth);
    dstruct.data.APwidth2 = cat(1,additionalData,dstruct.data.APwidth2);
    
    
    ostruct.handles(1).Children(end).Children(3).XData = ostruct.handles(1).Children(end).Children(3).XData - prop{1}(1).APt{end}(xth(1));
    ostruct.handles(1).Children(end).Children(2).XData = ostruct.handles(1).Children(end).Children(2).XData - prop{1}(2).APt{end}(xth(1));
    ostruct.handles(1).Children(end).Children(1).XData = ostruct.handles(1).Children(end).Children(1).XData - prop{1}(3).APt{end}(xth(1));
    
    ostruct.handles(1).Children(end-1).Children(3).XData = ostruct.handles(1).Children(end-1).Children(3).XData - prop{2}(1).APt{end}(xth(2));
    ostruct.handles(1).Children(end-1).Children(2).XData = ostruct.handles(1).Children(end-1).Children(2).XData - prop{2}(2).APt{end}(xth(2));
    ostruct.handles(1).Children(end-1).Children(1).XData = ostruct.handles(1).Children(end-1).Children(1).XData - prop{2}(3).APt{end}(xth(2));
    
    ostruct.handles(1).Children(end).YLim = [-80 70];
    ostruct.handles(1).Children(end).XLim = [-5 20];
    ostruct.handles(1).Children(end-1).YLim = [-80 70];
    ostruct.handles(1).Children(end-1).XLim = [-5 20];
    FontResizer
    % FigureResizer(5,8)
    % if ostruct.newborn
    %     tprint(fullfile(targetfolder_results,'RobustnessMatrix_newborn'),'-pdf');
    % else
    tprint(fullfile(targetfolder_results,sprintf('Fig.6_APwidthSpiking_%s',neuron.experiment)),'-pdf');
    
    % end
    ostruct.gap = 0.3;
    handle = Masterplotter(dstruct,gstruct,[],ostruct);
    figure(handle{1})
    tprint(fullfile(targetfolder_results,sprintf('Fig.6_APhalfwidth_%s',neuron.experiment)),'-pdf');
    figure(handle{2})
    tprint(fullfile(targetfolder_results,sprintf('Fig.6_APwidth_%s',neuron.experiment)),'-pdf');
    
    %% resonance test with Ba and ZD application, Figure_6-figure supplement2
    neuron = neuron_orig;
    neuron.experiment = strcat(neuron.experiment,'_resonance');
    holding_voltage = -95; % Stegen Hanuschkin Computer Sim 2012
    amp = 0.050; % 40 pA (less than +-50pA, Stegen 2012)
    neuron.params.tstop = 30000;
    neuron.params.celsius = 34.4;
    neuron.params.dt=0.5;
    
%     neuron = t2n_changemech(neuron,struct('gbar_HCN',2),'relative');
%     neuron = t2n_changemech(neuron,struct('gkbar_Kir21',2),'relative');
    figure;
    t2n_resonance(neuron,tree,amp,holding_voltage)
    
    neuron = t2n_blockchannel(neuron,'Kir21',100);
    t2n_resonance(neuron,tree,amp,holding_voltage)
    
    neuron = t2n_blockchannel(neuron,'HCN',100);
    t2n_resonance(neuron,tree,amp,holding_voltage)
    
    data{1} = importdata(fullfile(pwd,'raw data','StegenHanuschkin_Resonance_CTRL.csv'));
    data{2} = importdata(fullfile(pwd,'raw data','StegenHanuschkin_Resonance_Ba.csv'));
    data{3} = importdata(fullfile(pwd,'raw data','StegenHanuschkin_Resonance_Ba+ZD.csv'));
    plot(data{1}(:,1),data{1}(:,2),'k')
    plot(data{2}(:,1),data{2}(:,2),'b')
    plot(data{3}(:,1),data{3}(:,2),'r')
    ylim([0 350])
    
    FontResizer
    if isfield(ostruct,'figureheight') && isfield(ostruct,'figurewidth')
        FigureResizer(ostruct.figureheight,ostruct.figurewidth)
    end
        
    tprint(fullfile(targetfolder_results,sprintf('Fig6-Resonance_%s%s',neuron.experiment,str)),'-pdf');
    
    %%
end
%% here are simulations done in both mouse and rat
%% bAP simulation Figure 5, mimicking Krueppel et al 2011
neuron = neuron_orig;
ostruct.dist = 'Eucl.'; % PL., Eucl.
ostruct.relamp = 0;  % relative amplitudes
ostruct.plotData = ostruct.usemorph>=4; % plot experimental data
celsius_orig = neuron.params.celsius;
neuron.params.celsius = 33;  % temperature
neuron = t2n_Q10pas(neuron,neuron.params.celsius);
neuron.experiment = sprintf('%s_%d°',neuron.experiment,neuron.params.celsius);
if ~(ostruct.vmodel>=0)
    neuron = manipulate_Ra(neuron,1,'axon'); % necessary because otherwise the axon spikes permanently after the buzz and Ca decay measurement is not possible
    neuron.params.nseg = 1;  % necessary because dlambda calculation takes to long with high Ra
    neuron.params.accuracy = 1; % improve AIS segment number for more accurate simulation
    amp = 1700*0.001; %nA %AH99 model, cstep has to be bigger otherwise not all cells fire
else
    amp = 1300*0.001; %nA
end
t2n_bAP(neuron,tree,amp,targetfolder_data)  % do the simulation

neuron.params.celsius = celsius_orig; % back to normal temp
[bAPdisthm,mveloc_dend,mveloc_farax,mveloc_nearax,fig] = t2n_plotbAP(targetfolder_data,neuron,ostruct,targetfolder_results);
%%  Calcium dynamics Figure 5, mimicking Stocca et al 2008
neuron = neuron_orig;
ostruct.simple = 0;  % only recording along longest dendrite
ostruct.reduce = 0;  % only measure every third node
ostruct.dist = 'Eucl.'; % PL., Eucl.
ostruct.relamp = 0;  % relative amplitudes
celsius_orig = neuron.params.celsius;
neuron.params.celsius = 24;  % temperature
neuron = t2n_Q10pas(neuron,neuron.params.celsius);
neuron.experiment = sprintf('%s_%d°',neuron.experiment,neuron.params.celsius);
if ~(ostruct.vmodel>=0)
    ostruct.cai = 'caim_Caold';  % AH99 model, does not have an own calcium buffer mechanism, thus variable has a different name
    neuron = manipulate_Ra(neuron,1,'axon'); % necessary because otherwise the axon spikes permanently after the buzz and Ca decay measurement is not possible
    neuron.params.nseg = 1;  % necessary because dlambda calculation takes to long with high Ra
    neuron.params.accuracy = 1; % improve AIS segment number for more accurate simulation
    ostruct.cstep = 1700*0.001; %nA %AH99 model, cstep has to be bigger otherwise not all cells fire
else
    ostruct.cai = 'cai';
    ostruct.cstep = 1300*0.001; %nA
end
aGC_CaDyn(neuron,tree,targetfolder_data,ostruct)  % do the simulation

neuron.params.celsius = celsius_orig; % back to normal temp
aGC_plotCaDyn(targetfolder_data,targetfolder_results,neuron,ostruct);

%% Sensitivity Matrix Fig. 6 & Suppl. Fig.
changs = {'','cm','Ra','pas','Kir21','HCN','cAMP','na8st','Kv14','Kv21','Kv34','Kv42','Kv7','BK','SK','Cav12','Cav13','Cav22','Cav32','E_K','E_N_a','E_P_a_s','[Ca]o','temp'};
rmatrix = -Inf(14,numel(changs));

type = 'up'; % possible values: up down   ... regulate channels up or down
neuron = neuron_orig;
neuron.params.exchfolder = 't2nexchange_aGCmorphsim4';
neuron.experiment = strcat('robust_',neuron.experiment);

if strcmp(type,'up')
    amount = [-100,-10,20,4]; %some fixed change values for cm, ek, ena, eca, epas
else
    amount = [50,10,-20,1];
end
w = waitbar(0,'Matrix is calculating for some hours, please get yourself one or more snickers...');

for v = 1:numel(changs) %
    
    switch changs{v}
        case 'cm'
            neuron = t2n_blockchannel(neuron,'pas',amount(1),[],changs{v}); % specify pas channel
        case 'Ra'
            neuron = t2n_blockchannel(neuron,'pas',amount(1),[],changs{v}); % specify pas channel
        case 'E_K'
            for t = 1:numel(tree)
                neuron.mech{t}.all.k_ion.ek = neuron.mech{t}.all.k_ion.ek + amount(2);
            end
            neuron.experiment = strcat(neuron.experiment,'_ek');
        case 'E_N_a'
            for t = 1:numel(tree)
                neuron.mech{t}.all.na_ion.ena = neuron.mech{t}.all.na_ion.ena + amount(3);
            end
            neuron.experiment = strcat(neuron.experiment,'_ena');
        case 'E_P_a_s'
            for t = 1:numel(tree)
                fields = fieldnames(neuron.mech{t});
                for f = 1:numel(fields)
                    if isfield(neuron.mech{t}.(fields{f}),'pas')
                        neuron.mech{t}.(fields{f}).pas.e = neuron.mech{t}.(fields{f}).pas.e + amount(2);
                    end
                end
            end
            neuron.experiment = strcat(neuron.experiment,'_epas');
        case '[Ca]o'
            for t = 1:numel(tree)
                neuron.mech{t}.all.ca_ion.cao0 = amount(4);
            end
            neuron.experiment = strcat(neuron.experiment,'_eca');
        case 'temp'
            neuron.params.celsius = neuron.params.celsius + amount(2);
            neuron.experiment = strcat(neuron.experiment,'_temp');
        case 'cAMP'
            neuron = t2n_changemech(neuron,struct('cAMP_HCN',1),'absolute'); % 1µM cAMP
        otherwise
            if ~isempty(changs{v})
                neuron = t2n_blockchannel(neuron,changs{v},amount(1));
            end
            
    end
    %a
    ostruct.passtest = 'Std';
    ostruct.show = 0;
    [Rin, tau, ~, Vrest] = t2n_passTests(neuron,tree,targetfolder_results,ostruct);
    % caution, VoltageClamps add the series resistance of the electrode (10MOhm
    % by now). With passive, capacitance is higher..? But Rin is always same
    rmatrix(1,v) = mean(Vrest);
    rmatrix(2,v) = mean(Rin);
    rmatrix(3,v) = mean(tau);
    % b
    ostruct.holding_voltage = -80;
    ostruct.duration = 200;
    ostruct.amp = 90/1000;
    neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
    ostruct.coarse = 0.5;   % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
    ostruct.data = 2;
    rmatrix(5,v) = nanmean(t2n_findCurr(neuron,tree,'spike'));
    
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)
    props = t2n_APprop(targetfolder_data,neuron,0.09);
    rmatrix(4,v) = nanmean(cellfun(@(y) nanmean(y),props.APiv(1,:)));
    rmatrix(6,v) = nanmean(cellfun(@(y) nanmean(y),props.APamp(1,:)));
    rmatrix(7,v) = nanmean(cellfun(@(y) nanmean(y),props.APwidth(1,:)));
    rmatrix(8,v) = nanmean(cellfun(@(y) nanmean(y),props.fAHPabs(1,:)));
    % c
    ostruct.duration = 1000;%200;%200;
    ostruct.amp = [100,200]/1000; % nA
    t2n_currSteps(neuron,tree,targetfolder_data,ostruct)
    props = t2n_APprop(targetfolder_data,neuron,0.1);
    rmatrix(9,v) =  mean(cellfun(@(y) numel(y),props.APind(1,:))); %# APs@100
    rmatrix(11,v) =  nanmean(cellfun(@(y) mean(y),props.APISI(1,:))); % ISI mean
    isia = NaN(numel(tree),1);
    for t = 1:numel(tree)
        if numel(props.APISI{1,t}) > 1
            isia(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
            if isia(t) == 0
                isia(t) = NaN;
            end
        end
    end
    rmatrix(12,v) = nanmean(isia);
    
    %d
    props = t2n_APprop(targetfolder_data,neuron,0.2);
    rmatrix(10,v) = mean(cellfun(@(y) numel(y),props.APind(1,:)));%# APs@200
    %e
    ostruct.relamp = 0;
    t2n_bAP(neuron,tree,[],targetfolder_data,1)
    [bAPdisthm,mveloc_dend,mveloc_farax,mveloc_nearax] = t2n_plotbAP(targetfolder_data,neuron,ostruct,targetfolder_results);
    rmatrix(13,v) = nanmean(mveloc_dend);
    rmatrix(14,v) = nanmean(mveloc_farax);
    
    neuron = neuron_orig;
    neuron.experiment = strcat('robust_',neuron.experiment);
    neuron.params.exchfolder = 't2nexchange_aGCmorphsim4';
    waitbar(v/numel(changs),w)
end
close(w)
values = {'V_R_e_s_t','R_i_n','membr. time const.','voltage thresh. AP','current thresh. AP','AP amp','AP width','AP fAHP','#APs  100 pA','#APs  200 pA','ISI','adaptation ratio','dendr. veloc  mean','ax. veloc  far ax'};%,'bAP PL dist half-max'};

figure;
if any(imag(rmatrix(:)))
    [x,y] = find(imag(rmatrix));
    warning('Warning, found imaginary values for %s during change of %s',values{x},changs{y})
end
imagesc(real(rmatrix(:,2:end))./repmat(real(rmatrix(:,1)),1,size(rmatrix,2)-1)-1)
set(gca,'YTick',1:15)
set(gca,'YTickLabel',values)
set(gca,'XTick',1:numel(changs)-1)
set(gca,'XTickLabel',changs(2:end))

set(gca,'CLim',[-0.5 0.5])
set(gca,'Box','off')
colorbar
o.image = 0;
o.border = [8,2];
FigureResizer(8,20,[],o)
if ostruct.newborn
    save(fullfile(targetfolder_data,sprintf('Fig.6_SensitivityMatrix_%s_newborn_%s.mat',type,neuron.experiment)),'rmatrix')
    tprint(fullfile(targetfolder_results,sprintf('Fig.6_SensitivityMatrix_%s_newborn_%s',type,neuron.experiment)),'-pdf');
else
    save(fullfile(targetfolder_data,sprintf('Fig.6_SensitivityMatrix_%s_%s.mat',type,neuron.experiment)),'rmatrix')
    tprint(fullfile(targetfolder_results,sprintf('Fig.6_SensitivityMatrix_%s_%s',type,neuron.experiment)),'-pdf');
end
