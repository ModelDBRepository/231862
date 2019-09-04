function fig = aGC_instFIplot(targetfolder_data,neuron,ostruct)
%
if ~isfield(ostruct,'duration') && ostruct.dataset < 7
    ostruct.duration = 200;
end
if ~isfield(ostruct,'show')
    ostruct.show = 1:2;
end

load(t2n_catName(targetfolder_data,'Exp_Spiking',strcat(neuron.experiment,'.mat')),'voltVec','timeVec','numspikes','params','cstepsSpikingModel','tree','nneuron')

if any(ostruct.usemorph == [2,3,5,6])  % artificial cells
    modelcol = [0 1 0];
else
    modelcol = [0 0 1];
end
[exp_iclamp,cstepsSpiking,rate] = load_ephys(ostruct.dataset,'CClamp');

if numel(ostruct.amp)==1
    s2= 1;
    s1 = find(cstepsSpiking == ostruct.amp/1000);
else
    if isfield(ostruct,'ampprop')
        s2=find(cstepsSpikingModel == ostruct.ampprop/1000);
        s1=find(cstepsSpiking == ostruct.ampprop/1000);
    else
        s2=find(cstepsSpikingModel == 0.09);
        s1 = find(cstepsSpiking == 0.09);
    end
end
if isempty(exp_iclamp)
    s1 = [];
end
instFImodel = NaN(numel(tree),floor(ostruct.duration/100));

if ~isempty(s2)
    for t = 1:numel(tree)
        for w = 1:floor(ostruct.duration/100)
            instFImodel(t,w) = sum(diff(voltVec{t,s2} > 0,1,1) == -1 & timeVec{t,s2}(2:end)>=55+(w-1)*100 & timeVec{t,s2}(2:end)<55+(w)*100)/0.1;
        end
    end
    if ~isempty(s1)
        tvec = (1/rate:1/rate:size(exp_iclamp,1)/rate)';
        for t = 1:size(exp_iclamp,2)
            for w = 1:floor(ostruct.duration/100)
                instFIexp(t,w) = sum(diff(exp_iclamp(:,t,s1) > 0,1,1) == -1 & tvec(2:end)>=55+(w-1)*100 & tvec(2:end)<55+(w)*100)/0.1;
            end
        end
    end
    
    if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(1))
        fig(1) = ostruct.handles(1);
        figure(fig(1))
    else
        fig(1) = figure;hold all,
    end
    xvec = (1:floor(ostruct.duration/100)) * 100;
    if ~isempty(s1)
        hp = patch ([xvec (fliplr (xvec))], [(mean(instFIexp,1) + std(instFIexp,[],1)) (fliplr (mean(instFIexp,1) - std(instFIexp,[],1)))], [0 0 0]);
        hold on
        set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
        plot (xvec, mean(instFIexp,1), 'k')
    end
    errorbar(xvec,mean(instFImodel,1),std(instFImodel,[],1),'Color',modelcol)
    xlim([0 xvec(end)+100])
    set(gca,'XTick',[0,xvec])
    yl = get(gca,'YLim');
    ylim([0 yl(2)])
    xlabel('Intervals [ms]')
    ylabel('Instantaneous frequency [Hz]')
end