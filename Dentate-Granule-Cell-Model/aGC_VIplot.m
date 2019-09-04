function aGC_VIplot(targetfolder_data,experiment,options)

if nargin < 3
    options.dataset = 2;
end

if ~isfield(options,'show')
    options.show = 3;
end


figure;clf;hold all,
colo = colorme({'dim blue','dark red','light blue','red','turquois','orange'});
for o = 1:6
    [exp_vclamp,vsteps,rate] = load_ephys(o,'VClamp');

    tvec = 1/rate:1/rate:size(exp_vclamp,1)/rate;
        
    curr_mature = squeeze(mean(exp_vclamp(194*rate+1:204*rate+1,:,:),1));
    basl = squeeze(mean(exp_vclamp(94*rate+1:104*rate+1,:,:),1));
    
    curr_mature = curr_mature-repmat(curr_mature(:,find(vsteps>= -65,1,'first')),1,numel(vsteps));  % -65 war ca das holding pot
    
    spik_ind =cat(2,false(size(curr_mature,1),sum(vsteps<-80)),squeeze(any(exp_vclamp(tvec>0&tvec<204,:,vsteps>=-80) < -300,1))); %index which cells spike
    curr_mature(spik_ind) = NaN;
    
    if any(options.show == [1 3])
        
        mIV = nanmean(curr_mature,1);
        stdIV = nanstd (curr_mature,1);
        hp = patch ([(mIV + stdIV) (fliplr (mIV - stdIV))],[vsteps (fliplr (vsteps))], colo{o});
        set (hp, 'facealpha', 0.4, 'edgecolor', 'none')
        plot (mIV,vsteps, 'Color',colo{o},'LineWidth',3,'LineStyle','.')
    end
    
    xlabel('Measured Current [pA]')
end

for o = 1:6
    [exp_iclamp,csteps,rate] = load_ephys(o,'CClamp');
    tvec = 1/rate:1/rate:size(exp_iclamp,1)/rate;
    vamp = squeeze(mean(exp_iclamp(tvec>205&tvec<255,:,:),1)-0*mean(exp_iclamp(tvec<55,:,:),1));
    vamp(squeeze(any(exp_iclamp(tvec>55&tvec<255,:,:)>0,1))) = NaN;  % delete spikers
    plot(csteps,vamp','color',colo{o})
end

if any(options.show == [2 3])
    if exist(t2n_catName(targetfolder_data,'Exp_VoltSteps',experiment,'.mat'),'file')
        load(t2n_catName(targetfolder_data,'Exp_VoltSteps',experiment,'.mat'))
        
        line(zeros(1,numel(vstepsModel)),vstepsModel,'LineStyle','--','Color',[0.5 0.5 0.5])
        p = plot(steadyStateCurrVec-repmat(steadyStateCurrVec(find(vstepsModel>=-76,1,'first'),:),size(steadyStateCurrVec,1),1),vstepsModel);
        for t =1:size(steadyStateCurrVec,2)
            set(p(t),'color',tree{t}.col{1})
        end
        ylabel('Holding Voltage [mV] corrected')
    end
    if exist(t2n_catName(targetfolder_data,'Exp_Spiking',experiment,'.mat'),'file')
        load(t2n_catName(targetfolder_data,'Exp_Spiking',experiment,'.mat'))
        
        vamp = NaN(numel(cstepsSpikingModel),3);
        for s = 1:numel(cstepsSpikingModel)
            tvec = timeVec{1,s};
            exp_iclamp = cell2mat(voltVec(:,s)');
            if ~any(exp_iclamp(tvec>55&tvec<255,:)>0,1)
                vamp(s,:) = squeeze(mean(exp_iclamp(tvec>205&tvec<255,:),1)-0*mean(exp_iclamp(tvec<55,:),1));
            end
            
        end
        plot(cstepsSpikingModel*1000,vamp,'color',[1 0 1])
    end
end

title('dim blue: mat1, light blue mat2, turqois matBaCl, dark red: young1, light red young2, orange youngBaCl')
FontResizer
% FigureResizer(5,8)
% tprint(fullfile(targetfolder_results,strcat('Fig2-IV',neuron.experiment)),'-png')
% close(fig)
