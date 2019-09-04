function aGC_AHP_plot(targetfolder_data,neuron)

load(t2n_catName(targetfolder_data,'Exp_msAHP',neuron.experiment,'.mat'))

figure
tit = {'Control','XE991 (M Channel blocker)','Apamin (SK blocker)'};
for f=1:size(AHPvoltVec,1)
    
    for s = 1:size(AHPvoltVec,2)
        subplot(3,1,s)
        hold all
        title(tit{s})
        if ~isempty(AHPtimeVec{f,s})
            plot(AHPtimeVec{f,s},squeeze(AHPvoltVec{f,s}),'LineWidth',1.5,'Color',tree{f}.col{1})
        end
        line([320,350],[-80 -80],'LineStyle','--','Color','k','LineWidth',5)
        line([600,700],[-80 -80],'Color','k','LineWidth',5)
    end
    
end
linkaxes

for f=1:size(AHPvoltVec,1)
    for s = 1:size(AHPvoltVec,2)
        timewindow1 = AHPtimeVec{f,s} >= 320 & AHPtimeVec{f,s} <= 350;
        timewindow2 = AHPtimeVec{f,s} >= 600 & AHPtimeVec{f,s} <= 700;
        mAHP(s,f) = mean(AHPvoltVec{f,s}(timewindow1))-AHPvoltVec{f,s}(end);
        sAHP(s,f) = mean(AHPvoltVec{f,s}(timewindow2))-AHPvoltVec{f,s}(end);
    end
end

display(mAHP)
display(sAHP)