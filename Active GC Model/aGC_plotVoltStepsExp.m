function aGC_plotVoltStepsExp(loadingfile,targetfolder_results,ostruct)
if nargin < 3
    ostruct = [];
end

steps = -120-12.1;

if ~isfield(ostruct,'subtract_hv')
    ostruct.subtract_hv = 0;
end

load(loadingfile,'neuron')
[exp_vclamp,vsteps,rate] = load_ephys(ostruct.dataset,'VClamp');

if isfield(ostruct,'handles') && ishandle(ostruct.handles(1))
    fig(1) = ostruct.handles(1);
    figure(fig(1))
else
    fig(1) = figure; 
end
hold all
if isfield(ostruct,'handles') && numel(ostruct.handles)>1 && ishandle(ostruct.handles(2))
    fig(2) = ostruct.handles(2);
    figure(fig(2))
    ax=gca;
else
    fig(2) = figure;ax=axes;
end
hold all
if ostruct.subtract_hv
    basel = mean(exp_vclamp(0*rate+1:104*rate+1,:,:),1);
    exp_vclamp = exp_vclamp(:,:,:) - repmat(basel,size(exp_vclamp,1),1,1); % subtract current at baseline holding voltage (as Mongiat did)
end

str =  '';
figure(fig(1))
for f = 1:size(exp_vclamp,2)
    for s = 1:size(exp_vclamp,3)
        subplot(floor(sqrt(size(exp_vclamp,3))),ceil(sqrt(size(exp_vclamp,3))),s)
        hold all
        plot(1/rate:1/rate:size(exp_vclamp,1)/rate,squeeze(exp_vclamp(:,f,s)))
        
        ylabel('Current [pA]')
        xlabel('Time [ms]')
        ylim([-200,200])
        title(sprintf('VClamp % 4.4g mV%s',vsteps(s),str));
        xlim([0 300])
        if any(vsteps(s) == steps)
            plot(ax,1/rate:1/rate:size(exp_vclamp,1)/rate,squeeze(exp_vclamp(:,f,s)))
            ylabel('Current [pA]')
            xlabel('Time [ms]')
            ylim([-400,200])
            title(sprintf('VClamp % 4.4g mV%s',vsteps(s),str));
            xlim([0 300])
        end
    end
end

figure(fig(2))
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth)
if isfield(ostruct,'savename')
    if ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
    end
else
    tprint(fullfile(targetfolder_results,expcat('IV_dyn',neuron.experiment)),'-pdf');
end