function handles = aGC_synstim_plot(neuron,type,ostruct,targetfolder_data,targetfolder_results,showv,c)
if ostruct.newborn
    str = '_newborn';
    str2 = 'bliblablub';
else
   str = '';
   str2 = '_newborn_';
end
if isfield(ostruct,'handles')
    handles = ostruct.handles;
else
    handles = [];
end
if nargin < 8 || isempty(c)
    c = 2;
end
col = colorme(c+1);
folds = dir(targetfolder_data);
folds = {folds.name};
folds = folds(cellfun(@(x) isempty(strfind(x,str2)) & ~isempty(strfind(x,sprintf('Exp_%s_%s%s',type,neuron.experiment,str))),folds));
if isempty(folds)
    error('No experiment found')
else
    answer = listdlg('ListString',folds,'selectionmode','single','promptstring','Please select file to load','listsize',[300,200]);
end
if ~isempty(answer)
    data = load(fullfile(targetfolder_data,folds{answer}));
else
    return
end
% data.freq = [10,20,40,100,150]; % Hz later -> MHz
% dd0 = 0:10:100;  % µm
% dt0 = -50:10:50; % ms

warning('off','signal:findpeaks:largeMinPeakHeight')

counter = 1;
switch type
    case 'test'
        
    case 'TBS'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        for t = 1:numel(data.tree)-numel(indstim)
            p = plot(data.out.t,cat(2,data.out.record{t}.cell.v{data.recnode{t}}));
            set(p,'Color',data.tree{t}.col{1})
            p(1).LineStyle = ':';
            p(2).LineStyle = '-.';
            p(3).LineStyle = '-';
            
        end
        legend(p,'Soma','100µm','200µm')
        xlabel('Time [ms]')
        ylabel('Membrane potential [mV]')
        
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        hold all;
        for t = 1:numel(data.tree)-numel(indstim)
            p = plot(data.out.t,data.out.record{t}.Exp2Syn.i{data.thesesynids{t}(1)});
            set(p,'Color',data.tree{t}.col{1})
            
        end
        xlabel('Time [ms]')
        ylabel('Synaptic conductance [nA]')
        
        
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        hold all;
        for t = 1:numel(data.tree)-numel(indstim)
            p = plot(data.out.t,cat(2,data.out.record{t}.cell.cai{data.recnode{t}})*1000000);
            set(p,'Color',data.tree{t}.col{1})
        end
        xlabel('Time [ms]')
        ylabel('Calcium concentration [nM]')
        
    case 'white'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        for t = 1:numel(data.tree)
            hold all,
            subplot(3,1,1)
            plot(data.out.t,data.out.record{t}.cell.v{1})
        end
        subplot(3,1,2)
        f = ones(1,1/data.neuron{1}.params.dt)*data.neuron{1}.params.dt/1;
        m_spikeMat = filtfilt(f,1,sum(spikeMat,1));
        plot(tvec,m_spikeMat)
        subplot(3,1,3)
        t2n_plotRaster(spikeMat,tvec)
        
    case 'regular'
        %         handles(end+1) = figure; hold all,
        %         subplot(5,1,1:4)
        %         plot(data.out.t,data.out.record{t}.cell.v{1})
        %         subplot(5,1,5)
        % %         line(repmat(1/data.freq:1/data.freq:data.neuron{1}.params.tstop,2,1),repmat([0;1],1,data.neuron{1}.params.tstop*data.freq),'color','k','LineWidth',3)
        %         ylim([-0.5 1.5])
        %         set(gca,'ytick',[])
        %         xlim([0,data.neuron{1}.params.tstop])
        data.freqmodel = NaN(numel(data.freq),numel(data.tree)-1);
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter));
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        for f = 1:numel(data.freq)
            subplot(numel(data.freq),1,f)
            line(repmat(1000/data.freq(f):1000/data.freq(f):data.neuron{1}.params.tstop,2,1),repmat([0;1],1,data.neuron{1}.params.tstop*data.freq(f)/1000),'color','k','LineWidth',2)
            for t = 1:numel(data.tree)-1
                hold all
                %                 plot(data.out{f}.t,data.out{f}.record{t}.cell.v{1})
                [~,ind] = findpeaks(data.out{f}.record{t}.cell.v{1},'MinPeakHeight',0);
                data.freqmodel(f,t) = numel(ind)/data.neuron{1}.params.tstop*1000;
                line(repmat(data.out{f}.t(ind)',2,1),repmat([t;t+1],1,numel(ind)),'color','b','LineWidth',2)
                if data.freq(f) == 40
                    if numel(handles)>=counter && ishandle(handles(counter))
                        figure(handles(counter));
                    else
                        handles(counter) = figure;hold all;
                    end
                    line(repmat(1000/data.freq(f):1000/data.freq(f):data.neuron{1}.params.tstop,2,1),repmat([0;1],1,data.neuron{1}.params.tstop*data.freq(f)/1000),'color','k','LineWidth',2)
                    line(repmat(data.out{f}.t(ind)',2,1),repmat([t;t+1],1,numel(ind)),'color',col{c},'LineWidth',2)
                    xlim([0 500])
                    figure(handles(counter-1))
                end
            end

            xlim([0,data.neuron{1}.params.tstop])
            %             ylim([0,numel(data.tree)+1-1])
        end
        if any(data.freq == 40)
            figure(handles(counter))
            FontResizer
            FigureResizer(ostruct.figureheight,ostruct.figurewidth)
            counter = counter+1;
        end

        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        errorbar(data.freq,mean(data.freqmodel,2),std(data.freqmodel,[],2))
        line([0, data.freq(end)],[0, data.freq(end)],'LineStyle','--','Color',[0.5 0.5 0.5])
        ylim([0, data.freq(end)])
        xlim([0, data.freq(end)+5])
        xlabel('data.freq_i_n [Hz]')
        ylabel('data.freq_o_u_t [Hz]')
    case 'temporal'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        data.freqmodel = NaN(numel(data.freq),numel(data.dt0),numel(data.tree)-2);
        for f = 1:numel(data.freq)
            for n = 1:numel(data.dt0)
                ind = (f-1)*numel(data.dt0) + n;
                subplot(numel(data.freq),numel(data.dt0),ind)
                hold all
                line(repmat(1000/data.freq(f):1000/data.freq(f):data.neuron{1}.params.tstop,2,1),repmat([-1;0],1,floor(data.neuron{1}.params.tstop*data.freq(f)/1000)),'color','k','LineWidth',2)
                line(repmat(1000/data.freq(f):1000/data.freq(f):data.neuron{1}.params.tstop,2,1)+data.dt0(n),repmat([-2;-1],1,floor(data.neuron{1}.params.tstop*data.freq(f)/1000)),'color','k','LineWidth',2)
                for t = 1:numel(data.tree)-2
                    if showv
                        plot(data.out{ind}.t,data.out{ind}.record{t}.cell.v{1})
                        ylim([-80 numel(data.tree)])
                    else
                        ylim([0 numel(data.tree)]-2)
                    end
                    [~,ind2] = findpeaks(data.out{ind}.record{t}.cell.v{1},'MinPeakHeight',0);
                    data.freqmodel(f,n,t) = numel(ind2)/data.neuron{1}.params.tstop*1000;
                    line(repmat(data.out{ind}.t(ind2)',2,1),repmat([t-1;t],1,numel(ind2)),'color','b','LineWidth',2)
                    if data.freq(f) == 10
                        if numel(handles)>=counter && ishandle(handles(counter))
                            figure(handles(counter));
                        else
                            handles(counter) = figure;hold all;
                        end
                        if any(data.dt0(n) == [-15,0])
                            subplot(1,2,find(data.dt0(n)==[-15,0])), hold all
                            if t == 1
                                line(repmat(1000/data.freq(f):1000/data.freq(f):data.neuron{1}.params.tstop,2,1),repmat([-1;0],1,floor(data.neuron{1}.params.tstop*data.freq(f)/1000)),'color','k','LineWidth',2)
                                line(repmat(1000/data.freq(f):1000/data.freq(f):data.neuron{1}.params.tstop,2,1)+data.dt0(n),repmat([-2;-1],1,floor(data.neuron{1}.params.tstop*data.freq(f)/1000)),'color','k','LineWidth',2)
                            end
                            line(repmat(data.out{ind}.t(ind2)',2,1),repmat([t-1;t],1,numel(ind2)),'color',col{c},'LineWidth',2)
                            xlim([0 500])
                        end
                        figure(handles(counter-1))
                    end
                end
                xlim([0 data.neuron{1}.params.tstop])
            end
        end
        
        FontResizer
%         FigureResizer(ostruct.figureheight,ostruct.figurewidth,[],ostruct)
        if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
                tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
        else
            tprint(fullfile(targetfolder_results,sprintf('TemporalSumRaster_%s_%s%s',type,data.neuron{1}.experiment,str)),'-pdf');
        end
        if any(data.freq == 10)
            figure(handles(counter))
            FontResizer
%             FigureResizer(ostruct.figureheight,ostruct.figurewidth)
            counter = counter+1;
        end
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        col = colorme('dim blue','pink','cyan');
        for f = 1:numel(data.freq)
            subplot(numel(data.freq)/2,numel(data.freq)/2,f)
            hold all,
            title(sprintf('data.frequency %g Hz',data.freq(f)))
%             errorbar(data.dt0,mean(data.freqmodel(f,:,:),3)/data.freq(f),std(data.freqmodel(f,:,:),[],3)/data.freq(f),'color',col{f})
            plot(data.dt0,squeeze(mean(data.freqmodel(f,:,:),3))/data.freq(f),'color',col{ostruct.newborn+1})
            ylim([0 1])
            xlim([data.dt0(1) data.dt0(numel(data.dt0))])
            xlabel('\Deltat0 [ms]')
            ylabel('I/O data.freq. ratio')
        end
        FontResizer
%         FigureResizer(ostruct.figureheight,ostruct.figurewidth,[],ostruct)
        if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
                tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
        else
            tprint(fullfile(targetfolder_results,sprintf('TemporalSumPlot_%s_%s%s',type,data.neuron{1}.experiment,str)),'-pdf');
        end
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        imagesc(mean(data.freqmodel,3)./repmat(data.freq',1,numel(data.dt0)))
        set(gca,{'CLim','XTick','XTickLabel','YTick','YTickLabel','YDir'},{[0 1],1:numel(data.dt0),data.dt0,1:numel(data.freq),data.freq,'reverse'})
        xlabel('delta time [ms]')
        ylabel('input data.frequency [Hz]')
        colorbar
        FontResizer
%         ostruct.image = 1;
%         FigureResizer(5,6,[],ostruct);%ostruct.figureheight,ostruct.figurewidth,[],ostruct)
        if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
                tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
        else
            tprint(fullfile(targetfolder_results,sprintf('TemporalSumImag_%s_%s%s',type,data.neuron{1}.experiment,str)),'-pdf');
        end
    case 'spatial'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        data.freqmodel = NaN(numel(data.freq),numel(data.dd0),numel(data.tree)-1);

        for f = 1:numel(data.freq)
            for n = 1:numel(data.dd0)
                ind = (f-1)*numel(data.dd0) + n;
                subplot(numel(data.freq),numel(data.dd0),ind)
                hold all
                line(repmat(1000/data.freq(f):1000/data.freq(f):data.neuron{1}.params.tstop,2,1),repmat([-1;0],1,floor(data.neuron{1}.params.tstop*data.freq(f)/1000)),'color','k','LineWidth',2)
                for t = 1:numel(data.tree)-1
                    if showv
                    plot(data.out{ind}.t,data.out{ind}.record{t}.cell.v{1})
                    else
                        ylim([-1 numel(data.tree)])
                    end
                    [~,ind2] = findpeaks(data.out{ind}.record{t}.cell.v{1},'MinPeakHeight',0);
                    data.freqmodel(f,n,t) = numel(ind2)/data.neuron{1}.params.tstop*1000;
                    line(repmat(data.out{ind}.t(ind2)',2,1),repmat([t-1;t],1,numel(ind2)),'color','b','LineWidth',2)
                end
                xlim([0 data.neuron{1}.params.tstop])
            end
        end
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        hold all
        col = colorme(numel(data.freq));
        for f = 1:numel(data.freq)
            subplot(numel(data.freq),1,f)
%             errorbar(data.dt0,mean(data.freqmodel(f,:,:),3)/data.freq(f),std(data.freqmodel(f,:,:),[],3)/data.freq(f),'color',col{f})
            plot(data.dd0,squeeze(data.freqmodel(f,:,:))/data.freq(f),'color',col{f})
            ylim([0 1])
        end
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        imagesc(mean(data.freqmodel,3)./repmat(data.freq',1,numel(data.dd0)))
        set(gca,{'CLim','XTick','XTickLabel','YTick','YTickLabel'},{[0 1],1:numel(data.dd0),data.dd0,1:numel(data.freq),data.freq})
        xlabel('delta dist. [µm]')
        ylabel('input data.frequency [Hz]')
        colorbar
    case 'spatial2'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        for f = 1:numel(data.freq)
            for n = 1:numel(data.dd0)
                ind = (f-1)*numel(data.dd0) + n;
%                 ind2 = (f-1)*numel(data.dd0) + 1;
                subplot(numel(data.freq),numel(data.dd0),ind)
                hold all
%                 line(repmat(1000/data.freq(f):1000/data.freq(f):data.neuron{1}.params.tstop,2,1),repmat([-1;0],1,floor(data.neuron{1}.params.tstop*data.freq(f)/1000)),'color','k','LineWidth',2)
                for t = 1:numel(data.tree)-1
                    plot(data.out{ind}.t,data.out{ind}.record{t}.cell.v{1})
                end
                xlim([0 data.neuron{1}.params.tstop])
                ylim([-90 -50])
            end
        end
end
warning('on','signal:findpeaks:largeMinPeakHeight')

