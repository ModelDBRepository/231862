function aGC_AHP(neuron,tree,targetfolder_data)

LJP = 7;  % from the paper
cd(path)
neuron.params.accuracy = 1;  % for more nseg in axon and soma!
neuron.params.v_init = -58 - LJP;

neuron.params.tstop = 3000;
neuron.params.dt=0.05;
neuron.params.cvode = 1;
    
%!!!
meanhvol = -58 - LJP;   % corrected!!
neuron.params.skiprun = 0; %!!!!!!!!!
if ~exist('hstep2','var')
    hstep2 = t2n_findCurr(neuron,tree,meanhvol,[],'-q-d');
end
for t=1:numel(tree)
    neuron.APCount{t} = [1,-30];
    neuron.record{t}.cell = struct('record','v','node',1);
%     neuron.pp{t} = {'APCount',1,struct('thresh',-30)};
%     neuron.record{t} = cat(1,neuron.record{t},{1,'time','pp'});
end
nneuron = cell(3,1);
nneuron{1} = neuron;
for s = 1:3
    switch s
        case 1
            nneuron{s} = neuron;
            for t = 1:numel(tree)
                nneuron{s}.pp{t}.IClamp = struct('node',elecnode,'times',[-100,200,300],'amp', [hstep2(t) hstep2(t)+0.45 hstep2(t)]); %n,del,dur,amp
%                 nneuron{s}.tree = [1 2 3];
            end
        case 2
            nneuron{s} = nneuron{1};
            nneuron{s} = t2n_blockchannel(nneuron{s},'Kv7');
%             for t = 1:numel(tree)
%                 nneuron{s}.mech{t}.axonh.Kv723.gbar = nneuron{s}.mech{t}.axonh.Kv723.gbar / 10000;
%             end
        case 3
            nneuron{s} = nneuron{1};
            nneuron{s} = t2n_blockchannel(nneuron{s},'SK');
%             for t = 1:numel(tree)
%                 nneuron{s}.mech{t}.soma.SK2.gkbar = nneuron{s}.mech{t}.soma.SK2.gkbar /10000;
%                 nneuron{s}.mech{t}.GCL.SK2.gkbar = nneuron{s}.mech{t}.GCL.SK2.gkbar /10000;
%                 nneuron{s}.mech{t}.adendIML.SK2.gkbar = nneuron{s}.mech{t}.adendIML.SK2.gkbar /10000;
%                 nneuron{s}.mech{t}.adendOML.SK2.gkbar = nneuron{s}.mech{t}.adendOML.SK2.gkbar /10000;
%                 nneuron{s}.mech{t}.adendMML.SK2.gkbar = nneuron{s}.mech{t}.adendMML.SK2.gkbar /10000;
%             end
    end
    
end
nneuron = t2n_as(1,nneuron);

    [out, ~] = t2n(nneuron,tree,'-q-d-w');
    if isfield(out,'error')
        return
    end

for s = 1:3
    for t = 1:numel(tree)
        if isfield(out{s},'error')
            AHPvoltVec{t,s} = [] ;
            AHPtimeVec{t,s} = [];
        else
            AHPvoltVec{t,s} = out{s}.record{t}.cell.v{1} ;
            AHPtimeVec{t,s} = out{s}.t;
        end
    end
end

save(fullfile(targetfolder_data,'EphysModel',sprintf('Exp_msAHP_%s.mat',neuron.experiment)),'AHPvoltVec','AHPtimeVec','neuron','tree')
