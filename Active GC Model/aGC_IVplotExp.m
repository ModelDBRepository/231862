function fig = aGC_IVplotExp(ostruct,loadingfile)

load(loadingfile,'mholding_current','steadyStateCurrVec','holding_voltage','vstepsModel')
if ostruct.subtract_hv
    steadyStateCurrVec = steadyStateCurrVec - repmat(mholding_current,size(steadyStateCurrVec,1),1);
end
if nargin < 2
    ostruct.dataset = 2;
end
if ~isfield(ostruct,'single')
    ostruct.single = 0;
end
if ~isfield(ostruct,'subtract_hv')
    ostruct.subtract_hv = 0;
end
if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(1))
    fig(1) = ostruct.handles(1);
    figure(fig(1))
else
    fig(1) = figure;clf;hold all,
end

if ~all(ostruct.dataset == 0)
    [exp_vclamp,vsteps,rate] = load_ephys(ostruct.dataset,'VClamp');
    tvec = 1/rate:1/rate:size(exp_vclamp,1)/rate;
    if ~all(ostruct.dataset == 2.28)  % dont use that VClamp dataset, as it had been done after spiking experiment
        Rin2 = zeros(size(exp_vclamp,2),1);
        cap2 = Rin2;
        Rin = Rin2;
        cap = Rin2;
        for t = 1:size(exp_vclamp,2)
            d = exp_vclamp(:,t,vsteps==holding_voltage(1)+10);
            is = mean(d(tvec>190&tvec<204));
            i0 = mean(d(tvec<104));
            Rin2(t) = (10)/(is-i0)*1000;  % Rin (in MOhm) mV/pA
            x = tvec(d < is & tvec' > 100);
            y = d(d < is & tvec' > 100)-is;
            y = y(x<=x(1)+50);  % voltage step length of original capacitance measurement was only 50 ms hence cut vector thereafter
            x = x(x<=x(1)+50);
            cap2(t) = trapz(x,y)/(-10); % calculate capacitance as integral (charge) divided by amplitude of voltage step
            
            d = exp_vclamp(:,t,vsteps==holding_voltage(1)-10);
            is = mean(d(tvec>190&tvec<204));
            i0 = mean(d(tvec<104));
            Rin(t) = (-10)/(is-i0)*1000;  % Rin (in MOhm) mV/pA
            x = tvec(d < is & tvec' > 100);
            y = d(d < is & tvec' > 100)-is;
            y = y(x<=x(1)+50);  % voltage step length of original capacitance measurement was only 50 ms hence cut vector thereafter
            x = x(x<=x(1)+50);
            cap(t) = trapz(x,y)/(-10); % calculate capacitance as integral (charge) divided by amplitude of voltage step
        end
        fprintf('Mean Rin in Exp(@%gmV) is %g +- %g MOhm (s.e.m., -10mV)\n',holding_voltage(1)-10,mean(Rin),std(Rin)/sqrt(numel(Rin)))
        fprintf('Mean Rin in Exp(@%gmV) is %g +- %g MOhm (s.e.m., +10mV)\n',holding_voltage(1)+10,mean(Rin2),std(Rin2)/sqrt(numel(Rin2)))
        fprintf('\nMean capacitance in Exp(@%gmV) is %g +- %g pF (s.e.m. -10mV)',holding_voltage(1)-10,mean(cap),std(cap)/sqrt(numel(cap)))
        fprintf('\nMean capacitance in Exp(@%gmV) is %g +- %g pF (s.e.m. +10mV)\n',holding_voltage(1)+10,mean(cap2),std(cap2)/sqrt(numel(cap2)))
    end
    if exist('delind','var')
        meas_curr = squeeze(mean(exp_vclamp(194*rate+1:204*rate+1,setdiff(1:size(exp_vclamp,2),delind),:),1));
        basl = squeeze(mean(exp_vclamp(94*rate+1:104*rate+1,setdiff(1:size(exp_vclamp,2),delind),:),1));
    else
        meas_curr = squeeze(mean(exp_vclamp(194*rate+1:204*rate+1,:,:),1));
        basl = squeeze(mean(exp_vclamp(94*rate+1:104*rate+1,:,:),1));
    end
    
    if ostruct.subtract_hv
        meas_curr = meas_curr - basl;
    end
    
    gKirReal = zeros(size(meas_curr,1),1);
    Kirind_Mongiat = vsteps <= -122.1;
    otherind_Mongiat = vsteps >= -82.1 & vsteps <= -62.1;
    for t = 1:size(meas_curr,1)
        tmp = polyfit(vsteps(Kirind_Mongiat),meas_curr(t,Kirind_Mongiat),1) - polyfit(vsteps(otherind_Mongiat),meas_curr(t,otherind_Mongiat),1) ;
        gKirReal(t) = tmp(1);
    end
    mIV = mean(meas_curr,1);
    if ostruct.single
        plot(vsteps,meas_curr,'Color','k')
        vrest = zeros(size(meas_curr,1),1);
        for n = 1:size(meas_curr,1)
            vrest(n) = interp1(meas_curr(n,:),vsteps,0) ;
        end
        sprintf('Mean Vrest %g +- %g mV\n',mean(vrest),std(vrest))
    else
        
        stdIV = std (meas_curr,1);
        if ostruct.newborn
            if ostruct.dataset == 2.28 && exist(fullfile(pwd,'raw data','IV_29dpi_S3.csv'),'file')
                meas_curr = importdata(fullfile(pwd,'raw data','IV_29dpi_S3.csv'));
                mIV = meas_curr(1:19,2)';
                stdIV = meas_curr(20:end,2)' - mIV;
            else
                mIV = NaN(1,numel(vsteps));
                stdIV = NaN(1,numel(vsteps));
            end
        end
        hp = patch ([vsteps (fliplr (vsteps))], [(mIV + stdIV) (fliplr (mIV - stdIV))], [0 0 0]);
        set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
        if all(ostruct.show == 1)
            plot (vsteps, mIV,'Color',[0 0 0],'LineWidth',1,'LineStyle','-')
        end
    end
    fprintf('Kir Slope Conductance real cells %s\n',sprintf(' %.3g+-%.3g nS, ',mean(gKirReal),std(gKirReal)))
end
ylabel('Current [pA]')



%% add paper data
if ostruct.dataset==0
    mIV = mean(steadyStateCurrVec,2);
    
    Brenner05 = [-80,-20,600,10]; %HV, pA step, MOhm, stdevMOhm
    Mongiat09 = [-70-12.1,-10,224,7]; %HV, mV step, MOhm, stdevMOhm  LJP corrected
    SH07 = [-80,-3,308,26]; %HV, pA step, MOhm, stdevMOhm
    
    Mongiat09y = [-70-12.1,-10,519,30]; %HV, mV step, MOhm, stdevMOhm  LJP corrected
    %     Piatti11y21 = [-70-12.1,-10,636,94]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 21 dpi , LJP unknown
    Piatti11y28 = [-70-12.1,-10,667,67]; %HV, mV step, MOhm, stdevMOhm  LJP corrected  %28 dpi , LJP unknown
    %     Yang15y22a = [-70-12.1,-10,990,470]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 22 dpi Ascl , LJP unknown
    %     Yang15y22b = [-70-12.1,-10,1040,640]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 22 dpi Glast, LJP unknown
    %     Yang15y25a = [-70-12.1,-10,840,400]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 25 dpi Ascl , LJP unknown
    %     Yang15y25b = [-70-12.1,-10,790,450]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 25 dpi Glast, LJP unknown
    Yang15y28a = [-70-12.1,-10,560,200]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 28 dpi Ascl , LJP unknown
    Yang15y28b = [-70-12.1,-10,660,250]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 28 dpi Glast, LJP unknown
    
    
    SH04 = [-80,-5,232,78]; % lila %HV, mV step, MOhm, stdevMOhm  RAT!!!
    Staley92a = [-85,-15,107,NaN]; % rot %HV, mV step, MOhm, stdevMOhm  LJP corrected  % rat!
    Staley92b = [-85,+15,228,14.2]; %rot %HV, mV step, MOhm, stdevMOhm  LJP corrected % rat!
    MA14 = [-62,-50,230,15]; %orange %HV, pA step, MOhm, stdevMOhm  LJP corrected % rat!
    Mehranfard = [-70,-50,295.6,11.5]; %grün
    
    Brunner14yRAT21 = [-80,-10,378,20]; %HV, pA step, MOhm, stdevMOhm   % 21 dpi
    Brunner14yRAT28 = [-80,-10,358,13]; %HV, pA step, MOhm, stdevMOhm   % 26 dpi
    
    col = colorme({'red','yellow','dim green','violett','cyan','pink'});
    FigureResizer(ostruct.figureheight,ostruct.figurewidth)   % has to come before arrows to avoid distortions
    if ~ostruct.newborn
        if ostruct.usemorph <= 3  % mouse morphologies
            a=arrow([Brenner05(1),0+interp1(vstepsModel,mIV,Brenner05(1))],[Brenner05(1) + Brenner05(3) * Brenner05(2)/1000,Brenner05(2)+interp1(vstepsModel,mIV,Brenner05(1))]);
            set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            a=arrow([Mongiat09(1), interp1(vstepsModel,mIV,Mongiat09(1))],[Mongiat09(1) + Mongiat09(2), Mongiat09(2)/Mongiat09(3)*1000+interp1(vstepsModel,mIV,Mongiat09(1))]);
            set(a,'FaceColor',col{2},'EdgeColor',col{2},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            a=arrow([SH07(1),interp1(vstepsModel,mIV,SH07(1))],[SH07(1) + SH07(3) * SH07(2)/1000,SH07(2)+interp1(vstepsModel,mIV,SH07(1))]);
            set(a,'FaceColor',col{3},'EdgeColor',col{3},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        else                            % rat model
            a=arrow([Staley92a(1),interp1(vstepsModel,mIV,Staley92a(1))],[Staley92a(1) + Staley92a(2),Staley92a(2)/Staley92a(3)*1000+interp1(vstepsModel,mIV,Staley92a(1))]);
            set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            a=arrow([Staley92b(1),interp1(vstepsModel,mIV,Staley92b(1))],[Staley92b(1) + Staley92b(2),Staley92b(2)/Staley92b(3)*1000+interp1(vstepsModel,mIV,Staley92b(1))]);
            set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            a=arrow([MA14(1),interp1(vstepsModel,mIV,MA14(1))],[MA14(1) + MA14(3) * MA14(2)/1000 , MA14(2)+interp1(vstepsModel,mIV,MA14(1))]);
            set(a,'FaceColor',col{2},'EdgeColor',col{2},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            a=arrow([Mehranfard(1),interp1(vstepsModel,mIV,Mehranfard(1))],[Mehranfard(1) + Mehranfard(3) * Mehranfard(2)/1000 , Mehranfard(2)+interp1(vstepsModel,mIV,Mehranfard(1))]);
            set(a,'FaceColor',col{3},'EdgeColor',col{3},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            a=arrow([SH04(1), interp1(vstepsModel,mIV,SH04(1))],[SH04(1) + SH04(2), SH04(2)/SH04(3)*1000+interp1(vstepsModel,mIV,SH04(1))]);
            set(a,'FaceColor',col{4},'EdgeColor',col{4},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        end
    elseif ostruct.newborn
        if ostruct.usemorph <= 3  % mouse morphologies
            a=arrow([Mongiat09y(1), interp1(vstepsModel,mIV,Mongiat09y(1))],[Mongiat09y(1) + Mongiat09y(2), Mongiat09y(2)/Mongiat09y(3)*1000+interp1(vstepsModel,mIV,Mongiat09y(1))]);
            set(a,'FaceColor',col{2},'EdgeColor',col{2},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            
            %         a=arrow([Piatti11y21(1), interp1(vstepsModel,mIV,Piatti11y21(1))],[Piatti11y21(1) + Piatti11y21(2), Piatti11y21(2)/Piatti11y21(3)*1000+interp1(vstepsModel,mIV,Piatti11y21(1))]);
            %         set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            
            a=arrow([Piatti11y28(1), interp1(vstepsModel,mIV,Piatti11y28(1))],[Piatti11y28(1) + Piatti11y28(2), Piatti11y28(2)/Piatti11y28(3)*1000+interp1(vstepsModel,mIV,Piatti11y28(1))]);
            set(a,'FaceColor',col{3},'EdgeColor',col{3},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            %
            %         a=arrow([Yang15y22a(1), interp1(vstepsModel,mIV,Yang15y22a(1))],[Yang15y22a(1) + Yang15y22a(2), Yang15y22a(2)/Yang15y22a(3)*1000+interp1(vstepsModel,mIV,Yang15y22a(1))]);
            %         set(a,'FaceColor',col{4},'EdgeColor',col{4},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            %
            %         a=arrow([Yang15y22b(1), interp1(vstepsModel,mIV,Yang15y22b(1))],[Yang15y22b(1) + Yang15y22b(2), Yang15y22b(2)/Yang15y22b(3)*1000+interp1(vstepsModel,mIV,Yang15y22b(1))]);
            %         set(a,'FaceColor',col{4},'EdgeColor',col{4},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            %
            %         a=arrow([Yang15y25a(1), interp1(vstepsModel,mIV,Yang15y25a(1))],[Yang15y25a(1) + Yang15y25a(2), Yang15y25a(2)/Yang15y25a(3)*1000+interp1(vstepsModel,mIV,Yang15y25a(1))]);
            %         set(a,'FaceColor',col{5},'EdgeColor',col{5},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            %
            %         a=arrow([Yang15y25b(1), interp1(vstepsModel,mIV,Yang15y25b(1))],[Yang15y25b(1) + Yang15y25b(2), Yang15y25b(2)/Yang15y25b(3)*1000+interp1(vstepsModel,mIV,Yang15y25b(1))]);
            %         set(a,'FaceColor',col{5},'EdgeColor',col{5},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            %
            a=arrow([Yang15y28a(1), interp1(vstepsModel,mIV,Yang15y28a(1))],[Yang15y28a(1) + Yang15y28a(2), Yang15y28a(2)/Yang15y28a(3)*1000+interp1(vstepsModel,mIV,Yang15y28a(1))]);
            set(a,'FaceColor',col{6},'EdgeColor',col{6},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            
            a=arrow([Yang15y28b(1), interp1(vstepsModel,mIV,Yang15y28b(1))],[Yang15y28b(1) + Yang15y28b(2), Yang15y28b(2)/Yang15y28b(3)*1000+interp1(vstepsModel,mIV,Yang15y28b(1))]);
            set(a,'FaceColor',col{6},'EdgeColor',col{6},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            
        else
            a=arrow([Brunner14yRAT21(1),0+interp1(vstepsModel,mIV,Brunner14yRAT21(1))],[Brunner14yRAT21(1) + Brunner14yRAT21(3) * Brunner14yRAT21(2)/1000,Brunner14yRAT21(2)+interp1(vstepsModel,mIV,Brunner14yRAT21(1))]);
            set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
            
            a=arrow([Brunner14yRAT28(1),0+interp1(vstepsModel,mIV,Brunner14yRAT28(1))],[Brunner14yRAT28(1) + Brunner14yRAT28(3) * Brunner14yRAT28(2)/1000,Brunner14yRAT28(2)+interp1(vstepsModel,mIV,Brunner14yRAT28(1))]);
            set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        end
    end
end
