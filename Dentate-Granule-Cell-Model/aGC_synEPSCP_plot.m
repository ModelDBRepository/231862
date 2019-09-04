function aGC_synEPSCP_plot(targetfolder_data,targetfolder_results,neuron)




%    load(fullfile(targetfolder_data,'/Exp_synEPSCP_Ba.mat'))
%    sim1 = load(fullfile(targetfolder_data,'/Exp_synEPSCP_CTRL.mat'));
% else
    load(fullfile(targetfolder_data,sprintf('/Exp_synEPSCP_%s.mat',neuron.experiment)));
% end



figure;hold all;
%     p = plot(mean(abs(sim1.EPSC),2),mean(sim1.EPSP1,2));
%     errbar(p,cat(3,std(abs(sim1.EPSP1),[],2),std(sim1.EPSC,[],2)))
% end
p = plot(mean(abs(EPSC),2),mean(EPSP1,2));
errbar(p,cat(3,std(abs(EPSP1),[],2),std(EPSC,[],2)))
xlim([0 300]),ylim([0 35])
ylabel('peak EPSP [mV]')
xlabel('peak EPSC [pA]')

tprint(t2n_catName(targetfolder_results,'Fig.X-EPSPC',nneuron.experiment),'-HR-pdf');
