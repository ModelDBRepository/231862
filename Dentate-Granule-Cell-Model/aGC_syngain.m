function aGC_syngain(nneuron,tree,targetfolder_data,ostruct)

nneuron_orig = nneuron;
syn = 'Krueppel';
switch ostruct.mode
    case 1
        syn = 'alpha';
        tit = 'Alpha synapse CTRL';
    case 2
        tit = 'AMPAR+NMDAR CTRL';
    case 3
        tit = 'AMPAR (D-APV)';
    case 4
        tit = 'AMPAR+NMDAR (TTX)';
        nneuron = t2n_blockchannel(nneuron,'na8st',100);
    case 5
        tit = 'AMPAR+NMDAR (Nickel)';
        nneuron = t2n_blockchannel(nneuron,'Cav32',100);
end

rs = 15;
onset = 10;
nsyn = 10;

nneuron.params.v_init = -80;
nneuron.params.skiprun = 0;
% cstep = 1300*0.001; %nA
nneuron.params.tstop = 70;
% nneuron.params.dt=0.025;

ntree = numel(tree);
if strcmp(syn,'Krueppel')
    tree{end+1} = struct('artificial','NetStim','number',1,'start',10); % add a netstim to activate the krueppel synapses
    nneuron.params.cvode = 0;
    nneuron.params.dt = 0.05;
else
   nneuron.params.cvode = 1; 
end
for t=1:ntree
    ipar = ipar_tree(tree{t});
    ipar = ipar(T_tree(tree{t}),:);  % only paths from termination points
    ipar = ipar(tree{t}.R(ipar(:,1))~=find(strcmp(tree{t}.rnames,'axon')),:); % use only nonaxon terminals
    plen{t} = Pvec_tree(tree{t});
    eucl{t} = eucl_tree(tree{t});

    [~,ind] = min(abs(plen{t}(ipar(1,1:find(ipar(1,:)==0,1,'first')-1))-150));  % a point at dist 150
    thesesynids{t} = ipar(1,ind:ind+nsyn); % should result in 10 synapses with distance 1 µm (assuming internode distance of 1µm) at one branch
%     [tree{t},thesesynids{t}] = spines_tree (tree{t}, thesesynids{t}, 0.1, 0.5, 0.5, 0);
    
    switch syn
        case 'alpha'
            nneuron.pp{t}.AlphaSynapse = struct('node',thesesynids{t},'onset',onset,'tau',3,'gmax',0.0002);
        case 'Krueppel'
              nneuron.pp{t}.KrueppelAMPA = struct('node',thesesynids{t},'gmax',1.8e-5*25,'tau1',0.05,'tau2',2); % laut keller 1990 4 ms decay
              if ostruct.mode ~=3
                  nneuron.pp{t}.KrueppelNMDA = struct('node',thesesynids{t},'gmax',1.63e-5*25,'mg',2,'tau1',0.33,'tau2',50,'eta',0.05);%tau2b 254,..laut Keller ist rise time 4 to 9 ms..achtung das ist aber low-pass filtered
                  nneuron.con(t) = struct('source',struct('cell',ntree+1,'watch','on'),'target',struct('cell',t,'pp',{'KrueppelAMPA','KrueppelNMDA'},'node',thesesynids{t}),'delay',0,'threshold',0.5,'weight',1);
              else
                  nneuron.con(t) = struct('source',struct('cell',ntree+1,'watch','on'),'target',struct('cell',t,'pp','KrueppelAMPA','node',	thesesynids{t}),'delay',0,'threshold',0.5,'weight',1);
              end
    end
    nneuron.record{t}.cell = struct('node',[1,thesesynids{t}],'record',{'v','cai'});
end

hstep = t2n_findCurr(nneuron,tree,-82.1); %assuming a HP of xxx mV

% fig(1) = figure;hold all,
% fig(2) = figure;hold all,
for s = nsyn:-1:1
    neuron = nneuron;
    for t = 1:ntree
        switch syn
            case 'alpha'
                neuron.pp{t}.AlphaSynapse.node = neuron.pp{t}.AlphaSynapse.node(1:s);  % reduce number of synapses to current nsyns
            case 'Krueppel'
                neuron.pp{t}.KrueppelAMPA.node = neuron.pp{t}.KrueppelAMPA.node(1:s);  % reduce number of synapses to current nsyns
                if ostruct.mode ~=3
                    neuron.pp{t}.KrueppelNMDA.node = neuron.pp{t}.KrueppelNMDA.node(1:s);  % reduce number of synapses to current nsyns
                end
        end
        neuron.pp{t}.IClamp = struct('node',1,'times',-200,'amp',hstep(t)); %n,del,dur,amp
    end
    [out, ~] = t2n(neuron,tree,'-w-q-d');
    
%     figure(fig(1));p = plot(out.t,cat(2,out.record{t}.cell.v{thesesynids{t}}));
%     figure(fig(2));p = plot(out.t,cat(2,out.record{t}.cell.v{1}));
%     %     hold all;plot(out.t,out.record{t}.cell.v{thesesynids{t}(1)},'r')
% %     legend(p,sprintfc('%d µm',dists(1:numel(thesesynids{t}))))
%     ylabel('Membrane potential [mV]')
%     xlabel('Time [ms]')
    for t = 1:ntree
        EPSP(s,t) = max(out.record{t}.cell.v{1}(out.t< onset+20))-mean(out.record{t}.cell.v{1}(out.t<10));
    end
    
end
save(fullfile(targetfolder_data,sprintf('Exp_syngain_%s_%s.mat',tit,nneuron_orig.experiment)),'nneuron','tree','EPSP','nsyn','tit')