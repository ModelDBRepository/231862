function aGC_BKblockJaffe11(neuron,tree,targetfolder_results)

% cstep = 1.8; %nA !
cstep = 0.3; % 300 pA
neuron.params.accuracy = 1;  % for more nseg in axon and soma!
neuron.params.dt=0.05;
neuron.params.cvode = 1;
neuron.params.tstop = 400;    
neuron.params.skiprun = 0; %!!!!!!!!!
if ~exist('hstep','var')
    hstep = [];
end
hstep = t2n_findCurr(neuron,tree,-80,hstep,'-q-d');

for t = 1:numel(tree)
%     plen = Pvec_tree(tree{t});
%     nodes{t} = [find(plen == 0,1,'first'),find(plen >= 8 & tree{t}.R == find(strcmp(tree{t}.rnames,'axonh')),1,'first'),find(plen >= 100 & tree{t}.R ~= 1,1,'first')];
    nodes{t} = 1;
    neuron.record{t}.cell = struct('node',nodes{t},'record','v');
    neuron.APCount{t} = [1,-20];
%     neuron.stim{t} = {'IClamp',1,struct('times',[-400,30,30.8],'amp', [hstep(t) hstep(t)+cstep hstep(t)])}; %n,del,dur,amp
    neuron.pp{t}.IClamp = struct('node',1,'times',[-400,50,950],'amp', [hstep(t) hstep(t)+cstep hstep(t)]); %n,del,dur,amp
end
nneuron = cell(1,2);
for s = 1:2
    if s == 2
        nneuron{s} = t2n_blockchannel(neuron,'BK');
    else
        nneuron{s} = neuron;
    end
end
[out] = t2n(nneuron,tree,'-q-d-w');
% out = out{1};
if isfield(out,'error')
    return
end
col = colorme({'Black','Red'});
for t=1:numel(tree)
    figure,hold on
    for s = 1:2
        for f = 1:numel(nodes{t})
            plot(out{s}.t,out{s}.record{t}.cell.v{nodes{t}(f)},'Color',col{s})
        end
        xlim([29.5,35])
        try
            xlim([min(out{1}.APCtimes{t}{1}(10),out{2}.APCtimes{t}{1}(10))-1,max(out{1}.APCtimes{t}{1}(10),out{2}.APCtimes{t}{1}(10))+10]) % 10th spike!
        end
    end
    legend('Control','BK block')
    tprint(t2n_catName(targetfolder_results,sprintf('Fig5-blockBK%d',t),neuron.experiment),'-HR-pdf')
end