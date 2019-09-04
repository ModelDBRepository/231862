function aGC_APsomax(neuron,tree,treepath,targetfolder_results)

disp('echte current injection 2 oder weniger nA. und nur 0.5 ms duration\n')
disp('check das nochmal sobald passives modell steht. evt muss das axon weniger leaky sein?')
disp('check axon geschwindigkeit (teilweise von leakyness abh.) bzw axon AP broadness\n')

cstep = 1.8; %nA !
neuron.params.accuracy = 1;  % for more nseg in axon and soma!
neuron.params.dt=0.05;
neuron.params.cvode = 1;
neuron.params.tstop = 350;    
neuron.params.skiprun = 0; %!!!!!!!!!

for t = 1:numel(tree)
    tree{t} = bleb_tree(tree{t},30);
    ind = find(~cellfun(@isempty,strfind(tree{t}.rnames,'axon'))); % get axonal regions
    ind = find(any(repmat(tree{t}.R,1,numel(ind)) == repmat(ind,numel(tree{t}.R),1),2) & tree{t}.D == 2) ; % find the bleb
    plen = Pvec_tree(tree{t});
    nodes{t} = [1,ind,find(plen >= 100 & tree{t}.R ~= 1,1,'first')];
end
tree = t2n_writeTrees(tree,[],strcat(treepath(1:end-4),'_bleb.mtr'));
    
hstep = t2n_findCurr(neuron,tree,-80,[],'-q-d');

for t = 1:numel(tree)
%     eucl{t} = eucl_tree(tree{t});
    neuron.record{t}.cell = struct('node',nodes{t},'record','v');
    neuron.pp{t}.IClamp = struct('node',1,'times',[-400,30,30.8],'amp', [hstep(t) hstep(t)+cstep hstep(t)]); %n,del,dur,amp %prev 30.8
    
end


out = t2n(neuron,tree,'-q-d-w');

if isfield(out,'error') && out.error > 0
    return
end




col = colorme({'Black','Green','Red'});

for t=1%:numel(tree)
    figure,hold on
    for f = 1:numel(nodes{t})
        plot(out.t,out.record{t}.cell.v{nodes{t}(f)},'Color',col{f})
        
        vec = diff(out.record{t}.cell.v{nodes{t}(f)})./diff(out.t);
        ind = find(out.t > 30.8,1,'first');
        maxrise(f,t) = max(vec(ind:end));
        [~,ind] = max(out.record{t}.cell.v{nodes{t}(f)});
        maxdecay(f,t) = -min(vec(ind:end));
    end
    xlim([29.5,35])
    legend('Soma','Axon @ 30 µm','Dendrite 100µm')
    ylabel('Membrane voltage [mV]')
    xlabel('Time [ms]')
end
FontResizer
FigureResizer(5,8)
tprint(t2n_catName(targetfolder_results,'Fig.3-APinit',neuron.experiment),'-HR-pdf')


figure;
subplot(1,5,1)
p = errorbar([1,2],[297,mean(maxrise(1,:))],[12,std(maxrise(1,:))/sqrt(numel(tree))]);
set(p,'LineStyle','none','Color','k','Marker','o');
ylim([0 400]),xlim([0 3])
set(gca,{'XTick','XTickLabel'},{[1 2],{'Data','Model'}})
ylabel('Max. rate of rise(V s^-^1)')
title('Soma')

subplot(1,5,2)
p = errorbar([1,2],[485,mean(maxrise(2,:))],[12,std(maxrise(2,:))/sqrt(numel(tree))]);
set(p,'LineStyle','none','Color','r','Marker','o');
ylim([0 600]),xlim([0 3])
set(gca,{'XTick','XTickLabel'},{[1 2],{'Data','Model'}})
ylabel('Max. rate of rise(V s^-^1)')
title('Axon')

subplot(1,5,3)
p = errorbar([1,2],[70,mean(maxdecay(1,:))],[2,std(maxdecay(1,:))/sqrt(numel(tree))]);
set(p,'LineStyle','none','Color','k','Marker','o');
ylim([0 400]),xlim([0 3])
set(gca,{'XTick','XTickLabel'},{[1 2],{'Data','Model'}})
ylabel('Max. rate of decay(V s^-^1)')
title('Soma')

subplot(1,5,4)
p = errorbar([1,2],[59,mean(maxdecay(2,:))],[2,std(maxdecay(2,:))/sqrt(numel(tree))]);
set(p,'LineStyle','none','Color','r','Marker','o');
ylim([0 400]),xlim([0 3])
set(gca,{'XTick','XTickLabel'},{[1 2],{'Data','Model'}})
ylabel('Max. rate of decay(V s^-^1)')
title('Axon')

FontResizer
% FigureResizer(5,8)
tprint(t2n_catName(targetfolder_results,'Fig.3-APinitdata',neuron.experiment),'-HR-pdf')
