function maxdv = aGC_plotdVExp(targetfolder_results,ostruct)
checkthis = 90;
if ~isfield(ostruct,'savename')
    ostruct.savename = 'Test';
end
[exp_iclamp,cstepsSpiking,rate] = load_ephys(ostruct.dataset,'CClamp');
if ~isfield(ostruct,'amp')
    ostruct.amp = cstepsSpiking;
end

dfilt = diff(exp_iclamp,1,1);

ind = ismember(cstepsSpiking,ostruct.amp);
exp_iclampCut = exp_iclamp(:,:,ind);
cstepsSpikingCut = cstepsSpiking(ind);
dfiltCut = dfilt(:,:,ind);

if isfield(ostruct,'handles') && ishandle(ostruct.handles(1))
    fig(1) = ostruct.handles(1);
    figure(fig(1));ax(1) = gca; hold all
else
    fig(1) = figure('units','normalized','outerposition',[0 0 1 1]);ax(1) = axes;
end
if isfield(ostruct,'handles') && numel(ostruct.handles>1) && ishandle(ostruct.handles(2))
    fig(2) = ostruct.handles(2);
    figure(fig(1)), hold all
else
    fig(2) = figure;hold all
end
fig2(1) = figure;ax(2) = axes;hold all
xlabel('Cell Voltage [mV]')
fig2(2) = figure;ax(3) = axes;hold all
xlabel('Cell Voltage [mV]')



p = zeros(size(exp_iclampCut,2),1);
col = colorme(size(exp_iclampCut,2),'-grk');
figure(fig(1))
for s = 1:size(exp_iclampCut,3)
    subplot(floor(sqrt(size(exp_iclampCut,3))),ceil(sqrt(size(exp_iclampCut,3))),s)
    hold all
    for f=1:size(exp_iclampCut,2)
        this = exp_iclampCut(2:end,f,s);
        xlabel('Cell Voltage [mV] (corrected!)')
        ind = find(dfiltCut(:,f,s)>7,1,'first') + 5*rate;  % find first spike
        if isempty(ind)
            ind =1;
            linstyl = '-';
        else
            linstyl = ':';
        end
        dt = (diff(1/rate:1/rate:size(exp_iclampCut,1)/rate,1))';
        
        plot(this(ind:end),dfiltCut(ind:end,f,s)./dt(ind:end),'LineWidth',1.5,'Color',colorme(col{f},'brighter'),'LineStyle',linstyl);%exp_iclamp(:,f,s)))  % rest of spikes is black
        p(f) = plot(this(1:ind),dfiltCut(1:ind,f,s)./dt(1:ind),'LineWidth',1.5,'Color',col{f},'LineStyle','-');%exp_iclamp(:,f,s)))  % mark first spike red
        if cstepsSpikingCut(s) == checkthis/1000
            plot(ax(3),this(ind:end),dfiltCut(ind:end,f,s)./dt(ind:end),'LineWidth',1.5,'Color',col{f},'LineStyle',linstyl);%exp_iclamp(:,f,s)))  % rest of spikes is dashed
            plot(ax(2),this(1:ind),dfiltCut(1:ind,f,s)./dt(1:ind),'LineWidth',1.5,'Color',col{f},'LineStyle','-');%exp_iclamp(:,f,s)))  % first spike is straight line
        end
        ylabel('dV')
        xlim([-80 80])
        set(gca,'XTick',-80:40:80)
        ylim([-200 800])
        set(gca,'YTick',-200:200:800)
    end
    uistack((p),'top')
end

figure(fig(1))
xlim([-80 80])
set(gca,'XTick',-80:40:80)
ylim([-200 800])
set(gca,'YTick',-200:200:800)
tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotAll',ostruct.savename)),'-pdf')




figure(fig2(1))  % exp first spike
xlim([-80 80])
set(gca,'XTick',-80:40:80)
ylim([-200 800])
set(gca,'YTick',-200:200:800)
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotExp',ostruct.savename)),'-pdf')

figure(fig2(2))  %EXP 2nd spike
xlim([-80 80])
set(gca,'XTick',-80:40:80)
ylim([-200 800])
set(gca,'YTick',-200:200:800)
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotExp2',ostruct.savename)),'-pdf')

maxdv = NaN(size(exp_iclamp,2),size(exp_iclamp,3));
for s = 1:size(exp_iclamp,3)
    for f=1:size(exp_iclamp,2)
        maxdv(f,s) = max(dfilt(:,f,s))*rate;
    end
end
figure(fig(2))
hp = patch ([cstepsSpiking*1000 (fliplr (cstepsSpiking*1000))], [(mean(maxdv,1) + std(maxdv,[],1)) (fliplr (mean(maxdv,1) - std(maxdv,[],1)))], [0.6 0.6 0.6],'edgecolor','none');
p = plot (cstepsSpiking*1000, mean(maxdv,1), 'k');
uistack(p,'bottom')
uistack(hp,'bottom')
xlim([0 max(cstepsSpiking)*1000])
ylim([0 700])
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
ylabel('Maximal dV/dt [mV/ms]')
xlabel('current steps [pA]')
tprint(fullfile(targetfolder_results,sprintf('%s-MaxdV',ostruct.savename)),'-pdf')


fprintf('Max dv of experiment (@ %g pA): %g +- %g mV/ms (s.e.m.)\n',ostruct.ampprop*1000,mean(maxdv(:,ostruct.ampprop==cstepsSpikingCut)),std(maxdv(:,ostruct.ampprop==cstepsSpikingCut))/sqrt(size(maxdv,1)) )