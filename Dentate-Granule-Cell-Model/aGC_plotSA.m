function aGC_plotSA(loadfile,targetfolder_results)

load(loadfile)

fig(1) = figure;

for f=1:size(voltVec,1)
    
    hold all
    if ~isempty(timeVec{f})
        plot(timeVec{f},squeeze(voltVec{f}),'LineWidth',1.5,'Color',tree{f}.col{1})
        for w = 1:10
            freq(w,f) = sum(timespikes{f} >= 100*(w-1)+50 & timespikes{f} <= 100*w+50)/0.1; % calculate spiking frequency in 100 ms time windows
        end
    end
end
figure(fig(1))
xlabel('Time [ms]')
ylabel('Membrane potential [mV]')
FontResizer
FigureResizer(5,8)
tprint(t2n_catName(targetfolder_results,'Fig.4-SA',neuron.experiment),'-pdf')

fig(2) = figure;
hold on
for t = 1:numel(tree)
    plot(50:100:1000,freq(:,t),'Color',tree{t}.col{1},'Marker','x')
end
xlabel('Intervals [ms]')
ylabel('Frequency [Hz]')
FontResizer
FigureResizer(5,8)
tprint(t2n_catName(targetfolder_results,'Fig.4-SAfreq',neuron.experiment),'-pdf')
