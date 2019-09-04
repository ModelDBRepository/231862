function aGC_plotsAHP(targetfolder_data,targetfolder_results,neuron)

load(t2n_catName(targetfolder_data,'Exp_MA14stimsAHP',neuron.experiment,'.mat'))
col = colorme({'Black','Red'});
freq = zeros(2,numel(tree));
for t=1:numel(tree)
    figure,hold on
    for s = 1:2
        for f = 1:numel(nodes{t})
            [~,isi] = findpeaks(out{s}.record{t}.cell.v{nodes{t}(f)},'MinPeakHeight',10);
            isi = diff(out{s}.t(isi));
            freq(s,t) = mean(1./isi)*1000;
            plot(out{s}.t,out{s}.record{t}.cell.v{nodes{t}(f)},'Color',col{s},'linewidth',2)
        end
%         xlim([29.5,35])
    end
    xlabel('Time [ms]')
    ylabel('Membrane Potential [mV]')
    xlim([0 350])
    ylim([-90 50])
    FontResizer
    FigureResizer(5,8)
end
text(260,30,'Control','fontname','Arial','fontsize',8,'color',col{1})
text(260,15,sprintf('%.3g Hz',freq(1,t)),'fontname','Arial','fontsize',8,'color',col{1})
text(260,-5,'SK block','fontname','Arial','fontsize',8,'color',col{2})
text(260,-20,sprintf('%.3g Hz',freq(2,t)),'fontname','Arial','fontsize',8,'color',col{2})

tprint(fullfile(targetfolder_results,'Fig5-blockSK'),'-HR-pdf')