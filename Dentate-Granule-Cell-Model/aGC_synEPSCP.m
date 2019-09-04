function [EPSC,EPSP1] = aGC_synEPSCP(nneuron,tree,targetfolder_data)
dists = 0:50:300;
rs = 15;
onset = 10;
nsyns = 5:5:25;
nneuron.params.v_init = -80;
nneuron.params.skiprun = 0;
% cstep = 1300*0.001; %nA
nneuron.params.tstop = 50;
% neuron.params.dt=0.025;
nneuron.params.cvode = 1;

hstep = t2n_findCurr(nneuron,tree,-82.1); %assuming a HP of xxx mV


for t=1:numel(tree)
    ipar = ipar_tree(tree{t});
    ipar = ipar(T_tree(tree{t}),:);  % only paths from termination points
    plen{t} = Pvec_tree(tree{t});
    eucl{t} = eucl_tree(tree{t});
    %     ipar(ipar==0) = 1;
    % %     ipar = ipar(~any(tree{t}.R(ipar)==find(strcmp(tree{t}.rnames,'axon')),2),:); % delete paths which are axon
    %     if ostruct.simple
    %         uipar = ipar(1,:);
    %         nodes{t} = uipar(1:find(uipar==1,1,'first'));
    %     else
    %         nodes{t} = unique(ipar);
    %     end
    %     if ostruct.reduce
    %         nodes{t} = nodes{t}(1:3:end);
    %     end
    %     neuron.record{t}.cell = struct('node',nodes{t},'record','v');
    synids = find(tree{t}.R == find(strcmp(tree{t}.rnames,'adendMML')));  % all adendMML nodes
    thesesynids{t} = synids(randperm(numel(synids),nsyns(end)));
    ind = find(any(ipar == thesesynids{t}(1),2),1,'first')  ; % find terminal with a stimulated synapse on the way to soma
    dpath{t} = ipar(ind,1:find(ipar(ind,:)==0,1,'first')-1);
    [val,ia] = min(abs(repmat(plen{t}(dpath{t}),1,numel(dists)) - repmat(dists,numel(dpath{t}),1)),[],1);
    dpath{t} = dpath{t}(ia(val<5)); % get node indices at desired path lengths on that branch (max deviation from dist value = 5)
    
    nneuron.pp{t}.AlphaSynapse = struct('node',thesesynids{t},'onset',onset,'tau',3,'gmax',0.0002);
    
    nneuron.record{t}.cell = struct('node',[thesesynids{t}(1),dpath{t}],'record',{'v','cai'});
end
fig(1) = figure;hold all,
fig(2) = figure;hold all,

for s = numel(nsyns):-1:1
    neuron = nneuron;
    for t = 1:numel(tree)
        neuron.pp{t}.AlphaSynapse.node = neuron.pp{t}.AlphaSynapse.node(1:nsyns(s));  % reduce number of synapses to current nsyns       
        neuron.pp{t}.IClamp = struct('node',1,'times',-200,'amp',hstep(t)); %n,del,dur,amp
    end
    [out, ~] = t2n(neuron,tree,'-w-q-d');
    
    figure(fig(1));p = plot(out.t,cat(2,out.record{t}.cell.v{dpath{t}}));
    %     hold all;plot(out.t,out.record{t}.cell.v{thesesynids{t}(1)},'r')
    legend(p,sprintfc('%d µm',dists(1:numel(dpath{t}))))
    ylabel('Membrane potential [mV]')
    xlabel('Time [ms]')
    
    for t = 1:numel(tree)
        neuron.pp{t} = rmfield(neuron.pp{t},'IClamp');
        neuron.pp{t}.SEClamp = struct('node',1,'times',-200,'amp', -82.1,'rs',rs); %n,del,dur,amp
        neuron.record{t}.SEClamp = struct('node',1,'record','i');
    end
    [out2, ~] = t2n(neuron,tree,'-w-q-d');
    figure(fig(2));p = plot(out2.t,out2.record{t}.SEClamp.i{1}*1000);
    %     hold all;plot(out.t,out.record{t}.cell.v{thesesynids{t}(1)},'r')
    %     legend(p,sprintfc('%d µm',dists(1:numel(dpath{t}))))
    ylabel('Current [pA]')
    xlabel('Time [ms]')
    
    %     figure;p = plot(out2.t,cat(2,out2.record{t}.cell.v{dpath{t}}));
    %     legend(p,sprintfc('%d µm',dists(1:numel(dpath{t}))))
    %     ylabel('Membrane potential [mV]')
    %     xlabel('Time [ms]')
    for t = 1:numel(tree)
    EPSP1(s,t) = max(out.record{t}.cell.v{1}(out.t< onset+20))-mean(out.record{t}.cell.v{1}(out.t<10));
    EPSP2(s,t) = max(out2.record{t}.cell.v{1}(out2.t< onset+20))-mean(out2.record{t}.cell.v{1}(out2.t<10));
    EPSC(s,t) = 1000 * (min(out2.record{t}.SEClamp.i{1})-mean(out2.record{t}.SEClamp.i{1}(out2.t<10)));
    end
    
end

save(fullfile(targetfolder_data,sprintf('Exp_synEPSCP_%s.mat',neuron.experiment)),'nneuron','tree','EPSP1','EPSC','dists','nsyns')
