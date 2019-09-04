function aGC_syngain_plot(neuron,targetfolder_data,ostruct)

switch ostruct.mode
    case 1
        tit = 'Alpha synapse CTRL';
    case 2
        tit = 'AMPAR+NMDAR CTRL';
    case 3
        tit = 'AMPAR (D-APV)';
    case 4
        tit = 'AMPAR+NMDAR (TTX)';
    case 5
        tit = 'AMPAR+NMDAR (Nickel)';
end

load(fullfile(targetfolder_data,sprintf('Exp_syngain_%s_%s.mat',tit,neuron.experiment)))


figure;hold all
plot(EPSP(1)*(1:nsyn),EPSP,'Marker','x')
line([0 12],[0 12],'LineStyle','--','Color',[0.5 0.5 0.5])
ylabel('measured EPSP [mV]')
xlabel('arithmetic sum [mV]')
xlim([0 12])
ylim([0 20])
title(tit)