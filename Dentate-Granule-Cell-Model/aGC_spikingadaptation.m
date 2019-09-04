function aGC_spikingadaptation(neuron,tree,targetfolder_data,holding_voltage)
%
current = 150*0.001;


neuron.params.celsius = 33;
neuron.params.accuracy = 1;  % for more nseg in axon and soma!
neuron.params.tstop = 1500;
neuron.params.dt=0.05;
neuron.params.cvode = 1;

hstep = t2n_findCurr(neuron,tree,holding_voltage,[],'-q-d');

for t=1:numel(tree)
    neuron.APCount{t} = [1,-30];
end

nneuron{1} = neuron;
for t = 1:numel(tree)
    nneuron{1}.pp{t}.IClamp = struct('node',1,'times',[-100,50,1050],'amp', [hstep(t) hstep(t)+current hstep(t)]); %n,del,dur,amp
    nneuron{1}.record{t}.cell = struct('node',1,'record','v');
end


[out, ~] = t2n(nneuron,tree,'-q-d-w',exchfolder);
if out{1}.error
    return
end
voltVec = cell(numel(tree),1);
timeVec = voltVec;
timespikes = voltVec;
for t = 1:numel(tree)
    if isfield(out{1},'error') && out{1}.error > 0
        voltVec{t} = [] ;
        timeVec{t} = [];
        timespikes{t} = [];
    else
        voltVec{t,1} = out{1}.record{t}.cell.v{1} ;
        timeVec{t,1} = out{1}.t;
        timespikes{t} = out{1}.APCtimes{t}{1};
    end
end
save(fullfile(targetfolder_data,sprintf('Exp_Adaptation_%s.mat',neuron.experiment)),'voltVec','timeVec','timespikes','current','tree','neuron')
end