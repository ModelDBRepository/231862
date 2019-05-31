function fig = aGC_plotCaDyn(targetfolder_data,targetfolder_results,nneuron,ostruct)
style = {'-','--',':','-.'};
if ~iscell(nneuron)
    nneuron = {nneuron};
end

fig(1) = figure;hold on
fig(2) = figure;hold on
for n = 1:numel(nneuron)
    fig(n+3) = figure;clf,hold all
    axis off
end


for n = 1:numel(nneuron)
    load(t2n_catName(targetfolder_data,'Exp_CaDyn',nneuron{n}.experiment,'.mat'),'plotcaivals','nodes','neuron','tree','CaNodes','mCai','stdCai','mV','stdV','tim','tw','maxcai')
    
    spiked = cellfun(@(x) any(x>0),mV(:,1));
    
    
    
    
    xlims = [Inf Inf];
    ylims = xlims;
    for t = 1:numel(tree)
        rax = find(strncmp(tree{t}.rnames,'axon',4));
        dendind = find(all(repmat(tree{t}.R,1,numel(rax))~=repmat(rax,numel(tree{t}.R),1),2));
        axind = find(any(repmat(tree{t}.R,1,numel(rax))==repmat(rax,numel(tree{t}.R),1),2));
        if ostruct.show
            figure(fig(n+3));
            ptree = tran_tree(rot_tree(tran_tree(tree{t}),[],'-m3dY'),[350*t 300 0]);
            ptree.D(ptree.D<2) = 2;
            axind2 = intersect(nodes{t},axind);
            plotcaivals{t}(axind2) = NaN;
            xlims = [min(xlims(1),min(ptree.X(dendind))),max(xlims(2),max(ptree.X(dendind)))];
            ylims = [min(ylims(1),min(ptree.Y(dendind))),max(ylims(2),max(ptree.Y(dendind)))];
            plot_tree(ptree,plotcaivals{t});
        end
        
        for f =1:numel(CaNodes{t})
            if 1%t == 1
                figure(fig(2))
                hold on
                if ~all(stdCai{t,f}==0)
                    hp = patch ([tim', (fliplr (tim'))], [(mCai{t,f} + stdCai{t,f}), (fliplr (mCai{t,f} - stdCai{t,f}))],[0 0 0]);
                    
                    set (hp, 'facealpha', 0.2, 'edgecolor', 'none','facecolor',tree{t}.col{1});%col{f})
                end
                plot(tim,mCai{t,f},'Color',tree{t}.col{1},'LineStyle',style{f});%)
                figure(fig(1))
                hold on
                if ~all(stdV{t,f}==0)
                    hp = patch ([tim' (fliplr (tim'))], [(mV{t,f} + stdV{t,f}) (fliplr (mV{t,f} - stdV{t,f}))],[0 0 0]);
                    
                    set (hp, 'facealpha', 0.2, 'edgecolor', 'none','facecolor',tree{t}.col{1})
                end
                plot(tim,mV{t,f},'Color',tree{t}.col{1},'LineStyle',style{f})
                
            end
        end
    end
    figure(fig(n+3))
    ylim(ylims)
    xlim(xlims)
    ostruct.image = 1 ;
    FigureResizer(5,17,[],ostruct)
    c = colorbar;
%     c.Limits =[-80,80];
    %         pos = get(c,'Position');
    set(c,'Position',[0.93 0.35 0.02 0.4],'fontweight','bold','fontname','Arial')
    %         yl = get(c,'YLim');
%     set(c,'YTick',[-80,0,80])
    %         set(c,'YTick',[ceil(yl(1)),0,floor(yl(2))])
    tprint(t2n_catName(targetfolder_results,'CaMax-trees',nneuron{n}.experiment),'-SHR-tif')
    
    
    figure(fig(2))
    %     subplot(2,1,1)
    %     legend([p2{1}(end),p2{2}(end),p2{3}(end)],'Soma','Proximal','Distal')
    ylabel('Ca concentration [nM]')
    xlim([0 200])
    xlabel('Time [ms]')
    ylim([0, 1000])
    figure(fig(1))
    xlim([30 60])
    ylim([-80 50])
    ylabel('Membrane Potential [mV]')
    xlabel('Time [ms]')
    
    figure(fig(2))
    FontResizer
    FigureResizer(5,8)
    figure(fig(1))
    FontResizer
    FigureResizer(5,8)
    
    gstruct.ugroupdef = {{'Data Prox.','Proximal','Distal','Soma','Axon','Data MFB'}};
    if ~isempty(strfind(nneuron{n}.experiment,'_art'))
        gstruct.col = colorme('Black','Dim green','Dim green','Dim green','Dim green','Black');
    else
        gstruct.col = colorme('Black','Dim blue','Dim blue','Dim blue','Dim blue','Black');
    end
    gstruct.group{1} = 1:6;
    data_stocca = [0.23*1000,nanmean(nanmean(tw,3),1),43];
    sem = [0.03*1000,nanstd(nanmean(tw,3),[],1)/sqrt(numel(tree)),6];
    if ostruct.usemorph >= 4 && ~ostruct.newborn  % rat mature..add data
        data_stocca = data_stocca([1,4,5,3,2,6]);
        sem = sem([1,4,5,3,2,6]);
    else
        data_stocca = data_stocca([4,5,3,2]);
        sem = sem([4,5,3,2]);
        gstruct.ugroupdef = {gstruct.ugroupdef{1}(2:5)};
        gstruct.col = gstruct.col(2:5);
    end
    
    figure;%(fig(5))
    [~,~,barwidth] = barme([],data_stocca,sem,gstruct);
    ylabel('Ca^2^+ decay constant [ms]')
    ylim([0 300])
    FontResizer
    FigureResizer(5,8,[barwidth,0.4])
    tprint(t2n_catName(targetfolder_results,'CaDecay',nneuron{n}.experiment),'-pdf');
    
    data_stocca= [194,mean(nanmean(maxcai,3),1)*1e6,mean([0.91,1.16]*1000)];
    sem = [23,std(nanmean(maxcai,3),[],1)/sqrt(numel(tree))*1e6,std([0.91,1.16]*1000)/sqrt(2)];
    if ostruct.usemorph >= 4 && ~ostruct.newborn  % rat mature..add data
        data_stocca = data_stocca([1,4,5,3,2,6]);
        sem = sem([1,4,5,3,2,6]);
    else
        data_stocca = data_stocca([4,5,3,2]);
        sem = sem([4,5,3,2]);
    end
    
    figure;%(fig(2))
    [~,~,barwidth] = barme([],data_stocca,sem,gstruct);
    ylabel('Peak Ca amplitude [nM]')
    ax = gca;
    
    FontResizer
    if ax.YLim(2) > 2000
        ylim([0 3.5E+5])
    else
        ylim([0 1200])
    end
    FigureResizer(5,8,[barwidth,0.4])
    tprint(t2n_catName(targetfolder_results,'CaAmp',nneuron{n}.experiment),'-pdf');
    
    
    if ~all(spiked)
        display('CAUTION: Not all cells spiked!')
    end
    fprintf('Mean Calcium decay time: Axon: %g +- %g nM Soma: %g +- %g ms Proximal: %g +- %g ms Distal: %g +- %g ms\n',[mean(nanmean(tw,3),1);std(nanmean(tw,3),[],1)])
    fprintf('Mean Calcium peak amplitude was: Axon: %g +- %g nM Soma: %g +- %g nM Proximal: %g +- %g nM Distal: %g +- %g nM\n',[mean(nanmean(maxcai,3),1);std(nanmean(maxcai,3),[],1)]*1e6)
end