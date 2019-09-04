function fig = aGC_plotCurrStepsExp(targetfolder_data,targetfolder_results,neuron,ostruct,steps)

if nargin < 5
    steps = [0.03,0.075]; % 30 and 75 pA
end
if nargin < 4
    ostruct.dataset = 2;
end
load(t2n_catName(targetfolder_data,'Exp_Spiking',neuron.experiment,'.mat'))

[exp_iclamp,cstepsSpiking,rate] = load_ephys(ostruct.dataset,'CClamp');

if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(1))
    fig(1) = ostruct.handles(1);
    figure(fig(1))
    [~,ia] = intersect(round(cstepsSpiking*1000),round(cstepsSpikingModel*1000)); % to account for if the steps were not the same in both
    nnan = ~isnan(cstepsSpiking);
    if ~isempty(ia)
        exp_iclamp = exp_iclamp(:,:,ia);
        cstepsSpiking = cstepsSpiking(ia);
    end
else
    fig(1) = figure;hold all,
    nnan = ones(numel(cstepsSpiking),1);
end


for f = 1:size(exp_iclamp,2)
    for s = 1:size(exp_iclamp,3)
        if nnan(s)
            subplot(round(sqrt(size(exp_iclamp,3))),ceil(sqrt(size(exp_iclamp,3))),s)
            hold all
            if ostruct.dataset < 6
                ylim([-82 -42])
            end
            ylabel('Cell Voltage [mV] (corrected!)')
            title(sprintf('IClamp % 4.4g pA',cstepsSpiking(s)*1000));
            if ~all(isnan(exp_iclamp(:,f,s)))
                p = plot(1/rate:1/rate:size(exp_iclamp,1)/rate,squeeze(exp_iclamp(:,f,s)),'Color',[0.4 0.4 0.4]);%exp_iclamp_mature(:,f,s)))
                xlim([0 350])
                uistack(p,'bottom')
            end
        end
    end
end

if ~all(isnan(exp_iclamp))
    Rin = (max(exp_iclamp(:,:,cstepsSpiking==0.01),[],1)-mean(exp_iclamp(1:50*rate,:,cstepsSpiking==0.01),1))/0.01;
    fprintf('\nMean Rin in Exp(@+10pA) is %g +- %g MOhm (s.e.m., +10mV)\n',mean(Rin),std(Rin)/sqrt(numel(Rin)))
end

if exist('steps','var') && ~isempty(intersect(cstepsSpiking,steps))
    fig(2) = figure; hold all
    for ss = 1:numel(steps)
        s = find(cstepsSpiking == steps(ss));
        if ~isempty(s)
            if ostruct.newborn && ostruct.dataset == 2.28;
                ff = [2,7,12];
            else
                ff = [4 1 8]; 
            end
            for f=1:3%1:size(voltVec,1)
                p(f)= plot(1/rate:1/rate:size(exp_iclamp,1)/rate,squeeze(exp_iclamp(:,ff(f),s)),'LineWidth',1);%exp_iclamp_mature(:,f,s)))
            end
            set(p(1),'Color','k')
            set(p(2),'Color',[0.3 0.3 0.3])
            set(p(3),'Color',[0.6 0.6 0.6])
            uistack(p,'bottom')
        end
    end
    FontResizer
    FigureResizer(ostruct.figureheight,ostruct.figurewidth)
    xlabel('Time [ms]')
    ylabel('Membrane voltage [mV]')
    xlim([0 350])
    ylim([-85 -30])
    if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,strcat(ostruct.savename,'-SpikingExp')),'-pdf');
    else
        tprint(fullfile(targetfolder_results,strcat('Fig.2-SpikingExp_',neuron.experiment)),'-pdf');
    end
end