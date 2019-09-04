%% comments

%% initialize trees and mechanisms
%********* folders 
targetfolder_results = 'D:\GCModel\results';  % the folder where the graphs and pictures are saved
targetfolder_data = 'D:\GCModel\simdata';  % the folder where simulated data is saved (that you do not need to simulate it again if you just want to plot something)

ostruct = struct('plot','auto','show',3,'legend',0,'marker','o','sem',1,'FontSize',10,'FontType','Arial','barwidth',0.3,'figurewidth',8,'figureheight',5,'LineWidth',1,'grid',0,'priority','plot','ticklength',0.015);  % some options that are used when something is plotted
ostruct.usecol = 1;  % 0 = pseudorandom colors for each simulated cell in graph, 1 = green or blue grading
ostruct.vmodel = 1; % 0 = passive model, > 0 = active model, everything else (e.g. NaN) = old AH99 model
ostruct.scalespines = 1;  % scaling of g_pas and cm to implicitly model spines
ostruct.adjustloads = 0;  % the Hay et al 2013 implementation of adjust dendritic loads to reduce variability between cells (not used in publication)
ostruct.noise = 0;       % add noise to the membrane voltage by injecting gaussian noise current to the soma (not working with variable dt / cvode)
%
ostruct.reducecells = 0;  % reduce number of cells for faster simulation (e.g. for testing)
ostruct.usemorph = 1;  % 1 = all SH07, 2= synth mouseMat, 3= synth mouseYoung 4= Beining AAV, 5 = synth ratOld 6= synth ratYoung 7 = Claiborne,
ostruct.newborn = 0;  % 0 = adult GC model, 1 = young abGC model


if ostruct.usemorph >= 4
    ostruct.ratadjust = 1;  % adjust Kir channel in rats
else
    ostruct.ratadjust = 0;
end

%*****************************


if ostruct.newborn
%     % this is the good FI,deep fAHP one...
%     ostruct.channelblock = {'Kir21','Kv42','na8st','BK','Cav13','Kv21'};%,'Kv21','Kv42','Kv14','Kv34','na8st','Cav13'};%'Cav22'};     %{'Kv42','Kir21','Kv14','Kv34','pas','Kv21','na8st','Kv723'};%,'na8st','Kv723','na8st','Kv21'};%'except','Kir21'};%{'Kir','Kv42','HCN'};%'Kir','Kv42','Kv14','HCN'}; % Kir , SK etc
%     ostruct.blockamount = [73,50,25,60,50,75];%68 kir %,100,80,80,50,55,100];%  Kv42 50 Kv21 70
%      ostruct.specify = {'','','','','',''};

     %      % this is the used version
    ostruct.channelblock = {'Kir21','Kv42','na8st','BK','Cav13','Kv21','Kv723','BK'};%,'Cav22','Cav12','Cav32'};%,'Kv21','Kv42','Kv14','Kv34','na8st','Cav13'};%'Cav22'};     %{'Kv42','Kir21','Kv14','Kv34','pas','Kv21','na8st','Kv723'};%,'na8st','Kv723','na8st','Kv21'};%'except','Kir21'};%{'Kir','Kv42','HCN'};%'Kir','Kv42','Kv14','HCN'}; % Kir , SK etc
    ostruct.blockamount = [73,50,25,40,50,50,50,100];%,100,100,100];%68 kir %,100,80,80,50,55,100];%  Kv42 50 Kv21 70
     ostruct.specify = {'','','','gakbar','','','','gabkbar'};
    
    ostruct.scalespines = 0.3;  % means g_pas and cm are only scaled by 30% of the original spine densities due to reduced spine density in young abGCs
else
    ostruct.channelblock = {};
    ostruct.blockamount = [];
end


[tree,neuron,treeFilename] = GC_initModel(ostruct);  % initialize the model by loading the morphologies and setting the biophysical parameters

neuron_orig = neuron;

%% rewrite trees if necessary
origfolder = pwd;
cd(path)
tree = t2n_writeTrees(tree,[],treeFilename);
cd(origfolder)

%% plot trees
ostruct.show = 2;
ostruct.savename = 'Fig3_RatSynTrees';
t2n_plotTrees(tree,targetfolder_results,[],ostruct) % '-s'

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% experiments start here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% passive parameter tests (SH 2007)
neuron = neuron_orig;
ostruct.passtest = 'Mehranfard'; %Mongiat Mongiat2 SH Brenner
[Rin, tau, cap,Vrest] = t2n_passTests(neuron,tree,targetfolder_results,ostruct);
% caution, VoltageClamps add the series resistance of the electrode (10MOhm
% by now). With passive, capacitance is higher..? But Rin is always same

%% voltage clamp simulation (Mongiat 2009)
neuron = neuron_orig;

holding_voltage = -80; % mV
ostruct.subtract_hv = 1; % boolean subtract holding voltage current
ostruct.show = 1:2; % 0= nothing, 1 = only exp data, 2 = only model data
ostruct.coarse = 1;
neuron.params.cvode = 1;  % boolean if dt is constant (0) or variable (1)
neuron.params.dt = 0.25;  % this is ignored if cvode = 1
ostruct.single = 0;
%
t2n_voltSteps(neuron,tree,-120:10:-40,[],holding_voltage,targetfolder_data)
%
t2n_IVplot(targetfolder_data,neuron,ostruct);
if ostruct.dataset ~= 2.28
    t2n_plotVoltSteps(targetfolder_data,neuron,[],ostruct.subtract_hv);
end
%% current clamp simulation (Mongiat 2009, Brenner 2005)
neuron = neuron_orig;

ostruct.holding_voltage = NaN; %  % CAUTION! no LJP subtraction if concrete value is given
ostruct.duration = 200;%200;

ostruct.coarse = 0.5;

ostruct.amp = [20,65,80,120]/1000; %nA
neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
neuron.params.dt = 0.25;  % this is ignored if cvode = 1

ostruct.show = 1:2;
ostruct.ampprop = 65/1000;

%%%%%%%%%%%%%%%%%%%%%
% neuron.params.celsius = 33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t2n_currSteps(neuron,tree,targetfolder_data,ostruct);

% % VI plot
% % aGC_VIplot(targetfolder_data,neuron.experiment,ostruct)
% % Figure current clamp voltage traces
t2n_plotCurrSteps(targetfolder_data,targetfolder_results,neuron,ostruct);
% % FI curve
% 
% % if ostruct.maxamp >= 90
% % AP properties during current clamp
prop =     t2n_APprop(targetfolder_data,neuron,ostruct,tree);
% % end
t2n_FIplot(targetfolder_data,neuron,ostruct,targetfolder_results);

%% Compare spike soma/axon with a axon bleb (simulation!) SH2010
% aGC_APsomax(neuron,tree,fullfile(treepath,treename),targetfolder_results)



%% dV/V curves  ONLY WITH CONSTANT DT and with high time resolution!!!
neuron = neuron_orig;
neuron.experiment = strcat(neuron.experiment,'_dV');
ostruct.amp = [40,65,90,115]/1000; % nA
neuron.params.cvode = 0;  % boolean if dt is constant (0) or variable (1)
neuron.params.dt = 0.25;  % this is ignored if cvode = 1
ostruct.coarse = 0;
ostruct.show = 1:2;
ostruct.ampprop = 90/1000;
if ostruct.newborn
    ostruct.dataset =2.28;%4;
else
    ostruct.dataset =3;
end
t2n_currSteps(neuron,tree,targetfolder_data,ostruct)
t2n_plotdV(targetfolder_data,neuron,ostruct,targetfolder_results);

%% bAP simulation (Krueppel 2011)
neuron = neuron_orig;
% neuron.params.celsius = 33; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ostruct.dist = 'Eucl'; % PL, Eucl
ostruct.relamp = 0;

t2n_bAP(neuron,tree,[],targetfolder_data)

t2n_plotbAP(targetfolder_data,neuron,ostruct,targetfolder_results);

%% mAHP/sAHP simulation ********************************************************************

% aGC_AHP(neuron,tree,targetfolder_data)
%
% aGC_AHP_plot(targetfolder_data,neuron)

%% block BK (Brenner 2005, Jaffe 2011) -> broader AP
% aGC_BKblockJaffe11(neuron,tree,targetfolder_results)

%% get ratio of aBK to aBK+abBK (Shruti 2012)
% aGC_Shruti_current(neuron,tree,ostruct)
%!!!!!!!!!!!!!!!!!!!!!!

%% spiking adaptation, as in Mateos-Aparicio 2014 (SK,M-Current contribution) ACHTUNG RATTE !!!!!!!!!!
holding_voltage = -77; % mV
% ostruct.show = 1:2; % 1 = only exp data, 2 = only model data

aGC_spikingadaptation(neuron,tree,targetfolder_data,holding_voltage);
aGC_plotSA(t2n_catName(targetfolder_data,'Exp_Adaptation',neuron.experiment,'.mat'),targetfolder_results)
%% block SK (Mateos-Aparicio 2014) ACHTUNG RATTE !!!!!!!!!!
% aGC_sAHPstimMA14(neuron,tree,targetfolder_data)
% aGC_plotsAHP(targetfolder_data,targetfolder_results,neuron)

%% Calcium channel contributions (Eliot Johnston 1994)
aGC_Ca_proportions_Elliot(neuron,tree,targetfolder_results)

% %% dendritic synaptic stimulation EPSC/EPSP
% % 
% aGC_synEPSCP(neuron,tree,targetfolder_data,ostruct)
% aGC_synEPSCP_plot(targetfolder_data,targetfolder_results,neuron,ostruct);
% 
%% Krueppel linear gain
ostruct.mode = 3;  % 1 = Alpha syn, 2 = Krueppel AMPAR+NMDAR, 3 = D-APV NMDAR block, 4 = TTX, 5 = Nickel (T-type block)
aGC_syngain(neuron,tree,targetfolder_data,ostruct)
aGC_syngain_plot(neuron,targetfolder_data,ostruct)

%%
ostruct.handles = [];
type = 'temporal'; % test, theta-burst stimulation, regular input, temporal shift in input, spatial shift in input
ostruct.synmode = 1; % 1 = no adjustment of synapse number, 2 = newborns have only half of synapse number, 3 = adjust the number of synapses to get subthreshold response
% ostruct.figureheight = 3;
aGC_synstim(neuron,tree,type,ostruct,targetfolder_data)

ostruct.handles = aGC_synstim_plot(tree,type,ostruct,targetfolder_data,targetfolder_results,0);

%%
ostruct.newborn = 1;
aGC_synstim_plot(tree,type,ostruct,targetfolder_data,targetfolder_results,0,3);
aGC_synstim_plot(tree,type,ostruct,targetfolder_data,targetfolder_results,0,4);
ostruct.handles = []