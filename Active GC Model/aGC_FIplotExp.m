function fig = aGC_FIplotExp(targetfolder_data,targetfolder_results,neuron,ostruct)
%
if ~isfield(ostruct,'duration') && ostruct.dataset < 7
    ostruct.duration = 200;
end
if ~isfield(ostruct,'spikeThresh')
    ostruct.spikeThresh = -10; % std spike thresh of -10 mV
end

load(t2n_catName(targetfolder_data,'Exp_Spiking',neuron.experiment,'.mat'),'voltVec','timeVec','numspikes','cstepsSpikingModel','tree','nneuron')

[exp_iclamp,cstepsSpiking] = load_ephys(ostruct.dataset,'CClamp');

if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(2))
    fig(1) = ostruct.handles(1);
    figure(fig(1))
else
    fig(1) = figure;hold all,
end
if ostruct.duration == 200 && ostruct.dataset ~= 0
    x = cstepsSpiking*1000;
    mFI = mean (squeeze((sum(diff(exp_iclamp > ostruct.spikeThresh,1,1) == -1,1))));
    stdFI = std (squeeze((sum(diff(exp_iclamp > ostruct.spikeThresh,1,1) == -1,1))));
    mFI = mFI(1:numel(x));
    stdFI= stdFI(1:numel(x));
    hp = patch ([x (fliplr (x))], [(mFI + stdFI) (fliplr (mFI - stdFI))],[0.7 0.7 0.7],'edgecolor','none' );
    p = plot (x, mFI, 'k');
elseif ostruct.duration == 900 && exist(fullfile(pwd,'raw data','FI_Brenner.csv'),'file')
    dataBrenner = importdata(fullfile(pwd,'raw data','FI_Brenner.csv'));
    hp = patch ([50:50:400 (fliplr (50:50:400))], [(dataBrenner.data(2:2:end))' (fliplr (2*dataBrenner.data(1:2:end)'-dataBrenner.data(2:2:end)'))], [0 0 0]);
    set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
    p = plot (50:50:400, dataBrenner.data(1:2:end), 'k');
elseif ostruct.duration == 1000 && exist(fullfile(pwd,'raw data','FI_Mehranfard15b_RAT.csv'),'file')
    dataMA = importdata(fullfile(pwd,'raw data','FI_Mehranfard15b_RAT.csv'));
    MAstd = dataMA(7:end,2)-dataMA(1:6,2);
    hp(1) = patch ([dataMA(1:6,1); (flipud (dataMA(1:6,1)))], [dataMA(1:6,2)-MAstd ;(flipud (dataMA(1:6,2)+MAstd))], [0 0 0]);
    set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
    p(1) = plot (dataMA(1:6,1), dataMA(1:6,2), 'k');
    dataMA = importdata(fullfile(pwd,'raw data','FI_Mehranfard14_RAT.csv'));
    MAstd = -dataMA.data(1:2:end-1)+dataMA.data(2:2:end);
    hp(2) = patch ([(100:50:250)'; (fliplr (100:50:250)')], [dataMA.data(1:2:end-1)-MAstd ;(flipud (dataMA.data(1:2:end-1)+MAstd))], [0 0 0]);
    set (hp, 'facecolor',[0.7,0.7,0.7],'facealpha', 0.2, 'edgecolor', 'none')
    p(2) = plot (100:50:250, dataMA.data(1:2:end-1), 'Color',[0.7,0.7,0.7]);
    dataMA = importdata(fullfile(pwd,'raw data','FI_Mehranfard15a_RAT.csv'));
    MAstd = dataMA(2:2:end,2)-dataMA(1:2:end-1,2);
    hp(3) = patch ([(50:50:250)'; (flipud ((50:50:250)'))], [dataMA(1:2:end-1,2)-MAstd ;(flipud (dataMA(1:2:end-1,2)+MAstd))], [0 0 0]);
    set (hp, 'facecolor',[0.5,0.5,0.5],'facealpha', 0.2, 'edgecolor', 'none')
    p(3) = plot (50:50:250, dataMA(1:2:end-1,2),  'Color',[0.5,0.5,0.5]);
end
uistack(p,'bottom')
uistack(hp,'bottom')

if max(cstepsSpikingModel*1000) <= 120
    if ostruct.newborn
        ylim([0 10])
        set(gca,'YTick',0:2:10)
    else
        ylim([0 8])
        set(gca,'YTick',0:2:8)
    end
    xlim([0 120])
else
    xlim([40 310])
    ylim([0 60])
    set(gca,'YTick',0:20:60)
end
if isnan(ostruct.vmodel) && ostruct.usemorph >= 4 % AH99
    xlim([0 310])
    ylim([0 120])
    set(gca,'YTick',0:20:120)
end
xlabel('Current step [pA]')
ylabel('Number of spikes')
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth)
if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
    tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf')
else
    tprint(t2n_catName(targetfolder_results,'Fig.4-FI',neuron.experiment),'-pdf')
end



if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(2))
    fig(2) = ostruct.handles(2);
    figure(fig(2))
else
    fig(2) = figure;hold all,
end
if ostruct.duration == 200
    x = cstepsSpiking*1000;
    mFI = mFI / 200 * 1000; % number of spikes divided by time of stimulation (200 ms) = Hz
    stdFI = stdFI / ostruct.duration * 1000;
    hp = patch ([x (fliplr (x))], [(mFI + stdFI) (fliplr (mFI - stdFI))], [0.7 0.7 0.7],'edgecolor','none');
    hold on
    if all(ostruct.show == 1)
        p = plot (x, mFI, 'k');
    end
elseif ostruct.duration == 900
    hp = patch ([50:50:400 (fliplr (50:50:400))], [(dataBrenner.data(2:2:end))' (fliplr (2*dataBrenner.data(1:2:end)'-dataBrenner.data(2:2:end)'))]/0.9, [0 0 0]);
    set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
    p = plot (50:50:400, dataBrenner.data(1:2:end)/0.9, 'k');
elseif ostruct.duration == 1000 && exist(fullfile(pwd,'raw data','FI_Mehranfard15b_RAT.csv'),'file')
    dataMA = importdata(fullfile(pwd,'raw data','FI_Mehranfard15b_RAT.csv'));
    MAstd = dataMA(7:end,2)-dataMA(1:6,2);
    hp(1) = patch ([dataMA(1:6,1); (flipud (dataMA(1:6,1)))], [dataMA(1:6,2)-MAstd ;(flipud (dataMA(1:6,2)+MAstd))], [0 0 0]);
    set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
    p(1) = plot (dataMA(1:6,1), dataMA(1:6,2), 'k');
    dataMA = importdata(fullfile(pwd,'raw data','FI_Mehranfard14_RAT.csv'));
    MAstd = -dataMA.data(1:2:end-1)+dataMA.data(2:2:end);
    hp(2) = patch ([(100:50:250)'; (fliplr (100:50:250)')], [dataMA.data(1:2:end-1)-MAstd ;(flipud (dataMA.data(1:2:end-1)+MAstd))], [0 0 0]);
    set (hp, 'facecolor',[0.7,0.7,0.7],'facealpha', 0.2, 'edgecolor', 'none')
    p(2) = plot (100:50:250, dataMA.data(1:2:end-1), 'Color',[0.7,0.7,0.7]);
    dataMA = importdata(fullfile(pwd,'raw data','FI_Mehranfard15a_RAT.csv'));
    MAstd = dataMA(2:2:end,2)-dataMA(1:2:end-1,2);
    hp(3) = patch ([(50:50:250)'; (flipud ((50:50:250)'))], [dataMA(1:2:end-1,2)-MAstd ;(flipud (dataMA(1:2:end-1,2)+MAstd))], [0 0 0]);
    set (hp, 'facecolor',[0.5,0.5,0.5],'facealpha', 0.2, 'edgecolor', 'none')
    p(3) = plot (50:50:250, dataMA(1:2:end-1,2),  'Color',[0.5,0.5,0.5]);
end
uistack(p,'bottom')
uistack(hp,'bottom')

if ostruct.usemorph < 4 % mouse
    ylim([0 40])
    xlim([0 120])
else
    ylim([0 60])
    xlim([40 310])
end
if isnan(ostruct.vmodel) && ostruct.usemorph >= 4 % AH99
    xlim([0 310])
    ylim([0 120])
    set(gca,'YTick',0:20:120)
end
xlabel('current step [pA]')
ylabel('frequency [Hz]')
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth)
if isfield(ostruct,'savename')  && ~isempty(ostruct.savename)
    if ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,[ostruct.savename,'_Hz']),'-pdf')
    end
else
    tprint(t2n_catName(targetfolder_results,'Fig.4-FI_Hz',neuron.experiment),'-pdf')
end