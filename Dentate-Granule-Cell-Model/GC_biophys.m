function strct = GC_biophys(options)
% options:
% -p    include passive parameters (default)
% -a    include active parameters (default)
% -o    old AH99....also calcium channels with uncorrected Ca buffer shell model
% -s   orig na8st vshifts and rates as in SH10
% -n   latest changes as published in Beining et al 2017

%% Paramters
%Kir parameter
Kirusemodel = 'Own'; % 'Own', 'Hanuschkin'
%Kv4.2 parameter
Kv4model = {'Amarillo','Barghaan'};
Kv4composition = 4; % {'wt','wt + KChIP2','wt + DPPX','wt + KChIP2 + DPPX'};
Kv4usethismodel = 2;


Ca_tau = 240;  %calcium decay! orig 9 ms

if ~isempty(strfind(options,'-n'))
    on = struct('HCN',1,'Kir',1.07,'Kv11',0.5/4,'Kv14',0.5/4,'Kv21',0.6*1.3/1.1 ,'Kv33',0,'Kv34',5 /4/2*2.5,'Kv42',0.145 *3 *2/2/2 *2,'Kv72u3',6.7,'Na',1.728,'Nav19',0,'NCa',0.5,'LCa',6,'TCa',0.33,'BK',1*1.5*2*2/1.5,'SK',0.2/2);%*3
elseif ~isempty(strfind(options,'-k'))
    on = struct('HCN',1,'Kir',1.07*0.95,'Kv11',0.5/4,'Kv14',0.5/4,'Kv21',0.6*1.3 ,'Kv33',0,'Kv34',5 /4/2,'Kv42',0.145 *3 *2/2/2,'Kv72u3',6.7,'Na',1.728,'Nav19',0,'NCa',0.5,'LCa',6,'TCa',0.33,'BK',1*1.5*2*2,'SK',0.2/2);%*3
else
    on = struct('HCN',1,'Kir',1.07,'Kv11',0,'Kv14',0.5/4,'Kv21',0.6*1.3/1.1 ,'Kv33',0,'Kv34',5 /4/2,'Kv42',0.145 *3 *2/2/2,'Kv72u3',6.7,'Na',1.728,'Nav19',0,'NCa',0.5,'LCa',6,'TCa',0.33,'BK',1*1.5*2*2,'SK',0.2/2);
end

if nargin < 1 || isempty(options)
    options = '-a-p';
end

strct = struct;
if ~isempty(strfind(options,'-p'))
    if isempty(strfind(options,'-o')) % new model
        if ~isempty(strfind(options,'-y'))
            e_pas 	= -70;%-96.5  ; % chloride has ecl of -40mV !
            e_pas_axon = -70;
        else
            e_pas 	= -80;%-96.5  ; % chloride has ecl of -40mV ! TRP channels...
            e_pas_axon = -80;
        end
        % -80.4 for SH07 pure passive model reproduction
        Rm = 36.4;
        
        if ~isempty(strfind(options,'-n'))
            chg = 0.48;
            gax = chg*0.5;
            chg = chg * 1.05;
        else
            chg = 0.48;
            gax = chg;
        end
        if ~isempty(strfind(options,'-ra'))
            on.Kir = on.Kir * 2.5;2.7;%2.5
        end
        if ~isempty(strfind(options,'-cm'))
            cm = 0.7;
            Ra = 200; % 200 as SH2010 did...
        else
            cm = 0.9;
            Ra = 200; %200 as SH2010 did...
        end
        Raax = 100;
        
        if isempty(strfind(options,'-a'))  % do not change the passive parameters from their original values if not active model
            chg = 1;
            gax = 1;
        end
        strct.GCL.pas =        struct('cm',cm  ,'Ra', Ra,'g',(chg)/(Rm*1000),'e',e_pas); %SH07 cell4,6,7
        strct.adendIML.pas =   struct('cm',cm  ,'Ra', Ra,'g',(chg)/(Rm*1000),'e',e_pas); %SH07 cell4,6,7
        strct.adendMML.pas =   struct('cm',cm  ,'Ra', Ra,'g',(chg)/(Rm*1000),'e',e_pas); %SH07 cell4,6,7
        strct.adendOML.pas =   struct('cm',cm  ,'Ra', Ra,'g',(chg)/(Rm*1000),'e',e_pas); %SH07 cell4,6,7
        strct.soma.pas =       struct('cm',cm  ,'Ra', Ra,'g',chg/(Rm*1000),'e',e_pas); %SH07 cell4,6,7
        strct.axon.pas = struct('cm',cm, 'Ra', Raax, 'g',gax/(Rm*1000),'e',e_pas_axon);  % Ra of 120 because of SH10 saying in distal axon there is less Ra
        strct.axonh.pas = struct('cm',cm, 'Ra', Ra,  'g',gax/(Rm*1000),'e',e_pas_axon);
        
        
    else  % AH99 model
        e_pas = -70;
        strct.all.pas = struct('e',e_pas,'g',2.5e-5,'Ra',210,'cm',1);
        strct.adendIML.pas = struct('e',e_pas,'g',4e-5,'Ra',210,'cm',1.6);
        strct.adendMML.pas = struct('e',e_pas,'g',4e-5,'Ra',210,'cm',1.6);
        strct.adendOML.pas = struct('e',e_pas,'g',4e-5,'Ra',210,'cm',1.6);
    end
end



if ~isempty(strfind(options,'-a'))
    if isempty(strfind(options,'-o'))
        
        %%
        if on.Kir > 0
            Kirdistr = [1, 0.286, 0.143, 0.143, 0.143];
            % Kir in axons probably only for G-Protein activated Kirs (Kir3.x)...
            switch Kirusemodel
                case 'Own'
                    gkbar_Kir = 7e-4; %0.0007
                    % Kir current
                    spm_i = 3.5;%µM %7.5; %sollte über 3.5 sein für tau bei tiefen vclamps  %-10 possible!
                    chg_Kir = on.Kir * gkbar_Kir;%5.3; %2.5;
                    gsub = 0.15; % substate conductance (between 0.2 (spermidine) and 0.07 (spermine) Liu 2012)
                    fac = 0.001;
                    vshiftbs = 10;
                    vshiftbb = 0;
                    shiftmg = 1;
                    mg_i = 4;
                    b = 0.105;
                    fac = 0.005 ;
                    mg_i = 4;
                    if ~isempty(strfind(options,'-n'))
                        Kirdistr = [1.05, 1.05, 1.05, 1.05, 1.05, 0.5]*0.9;
                    else
                        Kirdistr = [1, 1, 1, 1, 1, 1]*0.9;
                    end
                    chg_Kir = chg_Kir *0.2;
                    
                    gsub = 0.25 ;
                    
                    spm_i = 1;
                    shiftmg = 0.5    ; % means only 0.5 of ek shift
                    vshiftbs = 0;
                    cas = 1/7;  % tau made larger due to kir 2.3
                    As = 0.2;
                    
                    strct.soma.Kir21 =     struct('gkbar',Kirdistr(1)*chg_Kir,'mg_i',mg_i,'spm_i',spm_i,'gsub',gsub,'fac',fac,'vshiftbs',vshiftbs,'vshiftbb',vshiftbb,'b',b,'shiftmg',shiftmg,'cas',cas,'As',As);
                    strct.GCL.Kir21 =      struct('gkbar',Kirdistr(2)*chg_Kir,'mg_i',mg_i,'spm_i',spm_i,'gsub',gsub,'fac',fac,'vshiftbs',vshiftbs,'vshiftbb',vshiftbb,'b',b,'shiftmg',shiftmg,'cas',cas,'As',As);
                    strct.adendIML.Kir21 = struct('gkbar',Kirdistr(3)*chg_Kir,'mg_i',mg_i,'spm_i',spm_i,'gsub',gsub,'fac',fac,'vshiftbs',vshiftbs,'vshiftbb',vshiftbb,'b',b,'shiftmg',shiftmg,'cas',cas,'As',As);
                    strct.adendMML.Kir21 = struct('gkbar',Kirdistr(4)*chg_Kir,'mg_i',mg_i,'spm_i',spm_i,'gsub',gsub,'fac',fac,'vshiftbs',vshiftbs,'vshiftbb',vshiftbb,'b',b,'shiftmg',shiftmg,'cas',cas,'As',As);
                    strct.adendOML.Kir21 = struct('gkbar',Kirdistr(5)*chg_Kir,'mg_i',mg_i,'spm_i',spm_i,'gsub',gsub,'fac',fac,'vshiftbs',vshiftbs,'vshiftbb',vshiftbb,'b',b,'shiftmg',shiftmg,'cas',cas,'As',As);
                    strct.axon.Kir21 = struct('gkbar',Kirdistr(end)*chg_Kir,'mg_i',mg_i,'spm_i',spm_i,'gsub',gsub,'fac',fac,'vshiftbs',vshiftbs,'vshiftbb',vshiftbb,'b',b,'shiftmg',shiftmg,'cas',cas,'As',As);
                    strct.axonh.Kir21 = struct('gkbar',Kirdistr(end)*chg_Kir,'mg_i',mg_i,'spm_i',spm_i,'gsub',gsub,'fac',fac,'vshiftbs',vshiftbs,'vshiftbb',vshiftbb,'b',b,'shiftmg',shiftmg,'cas',cas,'As',As);
                    
                case 'Hanuschkin'
                    Kirdistr(:) = 1;
                    vhalfl_kir_fit     = -107.609612;
                    kl_kir_fit         = 10;               %     // NOTE: sign changed in comparison with publication
                    vhalft_kir_fit     =   67.0828 ;
                    at_kir_fit         = 0.00610779;
                    bt_kir_fit         = 0.0817741  ;             %      // NOTE: error in publication (here the correct value)
                    gkbar_kir_fit = 9.1016e-5;  % this is what hanuschkin used in his paper
                    chg_Kir = on.Kir * gkbar_kir_fit;
                    strct.soma.HKir =     struct('gkbar',Kirdistr(1)*chg_Kir,'vhalfl',vhalfl_kir_fit,'kl',kl_kir_fit,'vhalft',vhalft_kir_fit,'at',at_kir_fit,'bt',bt_kir_fit);
                    strct.GCL.HKir =      struct('gkbar',Kirdistr(2)*chg_Kir,'vhalfl',vhalfl_kir_fit,'kl',kl_kir_fit,'vhalft',vhalft_kir_fit,'at',at_kir_fit,'bt',bt_kir_fit);
                    strct.adendIML.HKir = struct('gkbar',1.5*Kirdistr(3)*chg_Kir,'vhalfl',vhalfl_kir_fit,'kl',kl_kir_fit,'vhalft',vhalft_kir_fit,'at',at_kir_fit,'bt',bt_kir_fit);
                    strct.adendMML.HKir = struct('gkbar',2*Kirdistr(4)*chg_Kir,'vhalfl',vhalfl_kir_fit,'kl',kl_kir_fit,'vhalft',vhalft_kir_fit,'at',at_kir_fit,'bt',bt_kir_fit);
                    strct.adendOML.HKir = struct('gkbar',2.5*Kirdistr(5)*chg_Kir,'vhalfl',vhalfl_kir_fit,'kl',kl_kir_fit,'vhalft',vhalft_kir_fit,'at',at_kir_fit,'bt',bt_kir_fit);
            end
        end
        
        if on.HCN > 0
            e_hcn 		= -41.9;
            vhalfl_hcn 	= -100; % shifted hcn since is much more left in other models. Furthermore, I have much Mg2+ and I have implemented the cAMP dependency
            kl_hcn 	= 8;			%// NOTE: sign changed in comparison with publication
            vhalft_hcn 	=   30.4;			%// (see manuscript)
            at_hcn    	= 0.00052;
            bt_hcn     	= 0.2151;
            cAMP = 0; % mM cAMP concentration
            
            ghcnbar = 4e-6; %1.7561e-5;
            % according to Notomi und Shigemoto 2004 only in dendrites
            strct.adendIML.HCN = struct('cAMP',cAMP,'e',e_hcn,'gbar',on.HCN*ghcnbar,'vhalfl',vhalfl_hcn,'kl',kl_hcn,'vhalft',vhalft_hcn,'at',at_hcn,'bt',bt_hcn);
            strct.adendMML.HCN = struct('cAMP',cAMP,'e',e_hcn,'gbar',on.HCN*ghcnbar,'vhalfl',vhalfl_hcn,'kl',kl_hcn,'vhalft',vhalft_hcn,'at',at_hcn,'bt',bt_hcn);
            strct.adendOML.HCN = struct('cAMP',cAMP,'e',e_hcn,'gbar',on.HCN*ghcnbar,'vhalfl',vhalfl_hcn,'kl',kl_hcn,'vhalft',vhalft_hcn,'at',at_hcn,'bt',bt_hcn);
        end
        % Na channel
        if on.Na>0
            %**********
            if ~isempty(strfind(options,'-s'))
                addShift = 0;   %no additional shift in SH10
            else
                addShift = 10;  % additional shift in my model
            end
            if ~isempty(strfind(options,'-y'))
                vShift = 12+addShift ; %12 is default
                vShift_inact = 22-vShift;
            else
                vShift = 12+addShift ; %12 is default
                vShift_inact = 22-vShift;
            end
            %**********
            
            f = fopen('lib_mech/soma_st8.txt','r');
            rates_soma = textscan(f, '%.20f','Delimiter','\n');
            fclose(f);
            f = fopen('lib_mech/axon_st8.txt','r');
            rates_axon = textscan(f, '%.20f','Delimiter','\n');
            fclose(f);
            
            gnabar_soma = on.Na*0.03 * 1.7;  %conversion to siemens;
            gnabar_axonh= on.Na*0.2*1.5;
            gnabar_axon = on.Na*0.03 * 1.7;  %siemens;
            
            gnabar_dend = 0*on.Na*0.0005;  % only for test purpose
            if isempty(strfind(options,'-s'))
                rates_axon{1}(13:end) = [2.9807,0.4679,--0.0596,0.3962,2982.1,0.0635];%  ah = [0.3962,2982.1,0.0635]; bh = [2.9807,0.4679,--0.0596];
                rates_soma{1}(13:end) = [2.9713,0.6443,--0.0594,1.5860,2306.7,0.0493]; % new fitted inactivation kinetics since SH10 seems not to consider recov from inact values
            end
            strct.axon.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_axon{1}(1),'a1_1',rates_axon{1}(2),'b1_0',rates_axon{1}(3),'b1_1',rates_axon{1}(4),'a2_0',rates_axon{1}(5),'a2_1',rates_axon{1}(6),'b2_0',rates_axon{1}(7),'b2_1',rates_axon{1}(8),'a3_0',rates_axon{1}(9),'a3_1',rates_axon{1}(10),'b3_0',rates_axon{1}(11),'b3_1',rates_axon{1}(12),'bh_0',rates_axon{1}(13),'bh_1',rates_axon{1}(14),'bh_2',rates_axon{1}(15),'ah_0',rates_axon{1}(16),'ah_1',rates_axon{1}(17),'ah_2',rates_axon{1}(18),'gbar',gnabar_axon);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
            strct.axonh.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_axon{1}(1),'a1_1',rates_axon{1}(2),'b1_0',rates_axon{1}(3),'b1_1',rates_axon{1}(4),'a2_0',rates_axon{1}(5),'a2_1',rates_axon{1}(6),'b2_0',rates_axon{1}(7),'b2_1',rates_axon{1}(8),'a3_0',rates_axon{1}(9),'a3_1',rates_axon{1}(10),'b3_0',rates_axon{1}(11),'b3_1',rates_axon{1}(12),'bh_0',rates_axon{1}(13),'bh_1',rates_axon{1}(14),'bh_2',rates_axon{1}(15),'ah_0',rates_axon{1}(16),'ah_1',rates_axon{1}(17),'ah_2',rates_axon{1}(18),'gbar',gnabar_axonh);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
            strct.soma.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_soma);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
            
            if gnabar_dend > 0
                ax=1;  % give dendritic na8st the axon rates
                if ax
                    strct.GCL.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_axon{1}(1),'a1_1',rates_axon{1}(2),'b1_0',rates_axon{1}(3),'b1_1',rates_axon{1}(4),'a2_0',rates_axon{1}(5),'a2_1',rates_axon{1}(6),'b2_0',rates_axon{1}(7),'b2_1',rates_axon{1}(8),'a3_0',rates_axon{1}(9),'a3_1',rates_axon{1}(10),'b3_0',rates_axon{1}(11),'b3_1',rates_axon{1}(12),'bh_0',rates_axon{1}(13),'bh_1',rates_axon{1}(14),'bh_2',rates_axon{1}(15),'ah_0',rates_axon{1}(16),'ah_1',rates_axon{1}(17),'ah_2',rates_axon{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
                    strct.adendIML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_axon{1}(1),'a1_1',rates_axon{1}(2),'b1_0',rates_axon{1}(3),'b1_1',rates_axon{1}(4),'a2_0',rates_axon{1}(5),'a2_1',rates_axon{1}(6),'b2_0',rates_axon{1}(7),'b2_1',rates_axon{1}(8),'a3_0',rates_axon{1}(9),'a3_1',rates_axon{1}(10),'b3_0',rates_axon{1}(11),'b3_1',rates_axon{1}(12),'bh_0',rates_axon{1}(13),'bh_1',rates_axon{1}(14),'bh_2',rates_axon{1}(15),'ah_0',rates_axon{1}(16),'ah_1',rates_axon{1}(17),'ah_2',rates_axon{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
                    strct.adendMML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_axon{1}(1),'a1_1',rates_axon{1}(2),'b1_0',rates_axon{1}(3),'b1_1',rates_axon{1}(4),'a2_0',rates_axon{1}(5),'a2_1',rates_axon{1}(6),'b2_0',rates_axon{1}(7),'b2_1',rates_axon{1}(8),'a3_0',rates_axon{1}(9),'a3_1',rates_axon{1}(10),'b3_0',rates_axon{1}(11),'b3_1',rates_axon{1}(12),'bh_0',rates_axon{1}(13),'bh_1',rates_axon{1}(14),'bh_2',rates_axon{1}(15),'ah_0',rates_axon{1}(16),'ah_1',rates_axon{1}(17),'ah_2',rates_axon{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
                    strct.adendOML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_axon{1}(1),'a1_1',rates_axon{1}(2),'b1_0',rates_axon{1}(3),'b1_1',rates_axon{1}(4),'a2_0',rates_axon{1}(5),'a2_1',rates_axon{1}(6),'b2_0',rates_axon{1}(7),'b2_1',rates_axon{1}(8),'a3_0',rates_axon{1}(9),'a3_1',rates_axon{1}(10),'b3_0',rates_axon{1}(11),'b3_1',rates_axon{1}(12),'bh_0',rates_axon{1}(13),'bh_1',rates_axon{1}(14),'bh_2',rates_axon{1}(15),'ah_0',rates_axon{1}(16),'ah_1',rates_axon{1}(17),'ah_2',rates_axon{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
                else
                    strct.GCL.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
                    strct.adendIML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
                    strct.adendMML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
                    strct.adendOML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
                end
            end
%         else
%             % this is the old SH10 na8st channel distribution which is
%             % distance dependent and thus need the tree morphology for
%             % calculation of distribution. Is not used in the model
%             if isempty(tree)
%                 warning('Tree not given to GC_biophys. 8st sodium channel model cannot be inserted!')
%             else
%                 
%                 toaxon = false(numel(tree.X),1);
%                 toaxon(tree.R == find(strcmp(tree.rnames,'axon')) | tree.R == find(strcmp(tree.rnames,'axonh'))) = true;
%                 [tree2,order] = redirect_tree(tree,find(Pvec_tree(tree,toaxon)==1,1,'first'));  % redirect tree so that first axon point is root node (axon initial segment)
%                 plen = Pvec_tree(tree2); % get distances from axon initial segment
%                 plen(order) = plen; % get original order (from soma root node)
%                 
%                 gnabar = NaN(numel(tree.X),1);
%                 
%                 gnabar_axon= 38.6152 *nafac / 1000;  %conversion to siemens;
%                 gnabar_soma = 18.8 *nafac / 1000;  %conversion to siemens;
%                 gnabar_dend = 3.8 *nafac / 1000;%3.8 *nafac / 1000;
%                 a0 = 18.015027;
%                 gnabar(toaxon) = gnabar_soma + (gnabar_axon - gnabar_soma) .* (1-exp(-plen(toaxon)/5)) .* (1+ a0 * exp(-plen(toaxon)/10));
%                 gnabar(~toaxon) = gnabar_dend + (gnabar_soma - gnabar_dend) ./ (1+ exp((abs(plen(~toaxon))-80)/40));
%                 gnabar(~toaxon & gnabar == gnabar_dend) = NaN;
%                 
%                 gnabar(toaxon & gnabar == gnabar_axon) = NaN;
%                 gnabar(~toaxon & gnabar == gnabar_soma) = NaN;%
%                 
%                 rates = NaN(numel(rates_axon{1}),numel(tree.X));
%                 for r = 1:numel(rates_axon{1})
%                     rates(r,toaxon) = rates_axon{1}(r) - (rates_axon{1}(r) - rates_soma{1}(r))./(1+exp((plen(toaxon)-2)/2));
%                     rates(r,~toaxon) = rates_axon{1}(r) - (rates_axon{1}(r) - rates_soma{1}(r))./(1+exp((-plen(~toaxon)-2)/2));
%                     % ranged rates outside soma can be
%                     % ignored anyways because lambda of
%                     % change is only 2 µm --> rate soma is
%                     % reached already in the soma
%                     if any(strcmp(tree.rnames,'axonh'))
%                         rates(r,~(tree.R == find(strcmp(tree.rnames,'axonh'))) & ~(tree.R == find(strcmp(tree.rnames,'axon'))) & ~(tree.R == find(strcmp(tree.rnames,'soma')))) = NaN;  % really...no Na in dendrites
%                     else
%                         rates(r,~(tree.R == find(strcmp(tree.rnames,'axon'))) & ~(tree.R == find(strcmp(tree.rnames,'soma')))) = NaN;  % really...no Na in dendrites
%                     end
%                 end
%                 strct.range.na8st = struct('a1_0',rates(1,:),'a1_1',rates(2,:),'b1_0',rates(3,:),'b1_1',rates(4,:),'a2_0',rates(5,:),'a2_1',rates(6,:),'b2_0',rates(7,:),'b2_1',rates(8,:),'a3_0',rates(9,:),'a3_1',rates(10,:),'b3_0',rates(11,:),'b3_1',rates(12,:),'bh_0',rates(13,:),'bh_1',rates(14,:),'bh_2',rates(15,:),'ah_0',rates(16,:),'ah_1',rates(17,:),'ah_2',rates(18,:),'gbar',gnabar);
%                 strct.axon.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_axon{1}(1),'a1_1',rates_axon{1}(2),'b1_0',rates_axon{1}(3),'b1_1',rates_axon{1}(4),'a2_0',rates_axon{1}(5),'a2_1',rates_axon{1}(6),'b2_0',rates_axon{1}(7),'b2_1',rates_axon{1}(8),'a3_0',rates_axon{1}(9),'a3_1',rates_axon{1}(10),'b3_0',rates_axon{1}(11),'b3_1',rates_axon{1}(12),'bh_0',rates_axon{1}(13),'bh_1',rates_axon{1}(14),'bh_2',rates_axon{1}(15),'ah_0',rates_axon{1}(16),'ah_1',rates_axon{1}(17),'ah_2',rates_axon{1}(18),'gbar',gnabar_axon);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
%                 strct.axonh.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_axon{1}(1),'a1_1',rates_axon{1}(2),'b1_0',rates_axon{1}(3),'b1_1',rates_axon{1}(4),'a2_0',rates_axon{1}(5),'a2_1',rates_axon{1}(6),'b2_0',rates_axon{1}(7),'b2_1',rates_axon{1}(8),'a3_0',rates_axon{1}(9),'a3_1',rates_axon{1}(10),'b3_0',rates_axon{1}(11),'b3_1',rates_axon{1}(12),'bh_0',rates_axon{1}(13),'bh_1',rates_axon{1}(14),'bh_2',rates_axon{1}(15),'ah_0',rates_axon{1}(16),'ah_1',rates_axon{1}(17),'ah_2',rates_axon{1}(18),'gbar',gnabar_axon);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
%                 strct.soma.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_soma);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
%                 strct.GCL.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
%                 strct.adendIML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
%                 strct.adendMML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
%                 strct.adendOML.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
%                 strct.adendOMLout.na8st = struct('vShift',vShift,'vShift_inact',vShift_inact,'a1_0',rates_soma{1}(1),'a1_1',rates_soma{1}(2),'b1_0',rates_soma{1}(3),'b1_1',rates_soma{1}(4),'a2_0',rates_soma{1}(5),'a2_1',rates_soma{1}(6),'b2_0',rates_soma{1}(7),'b2_1',rates_soma{1}(8),'a3_0',rates_soma{1}(9),'a3_1',rates_soma{1}(10),'b3_0',rates_soma{1}(11),'b3_1',rates_soma{1}(12),'bh_0',rates_soma{1}(13),'bh_1',rates_soma{1}(14),'bh_2',rates_soma{1}(15),'ah_0',rates_soma{1}(16),'ah_1',rates_soma{1}(17),'ah_2',rates_soma{1}(18),'gbar',gnabar_dend);%,'maxrate',8000,'vshift_inact',vshift_inact,'vShift',vdonnan,);
%             end    
        end
        if on.Nav19 > 0   % persistent sodium channel (as Krueppel 2011 proposed)
            vshift = -18;
            strct.GCL.Nav19 = struct('gbar',  on.Nav19 * 0.0002,'vshift',vshift);
            strct.adendIML.Nav19 = struct('gbar',  on.Nav19 * 0.0002,'vshift',vshift);
            strct.adendMML.Nav19 = struct('gbar',  on.Nav19 * 0.0002,'vshift',vshift);
            strct.adendOML.Nav19 = struct('gbar',  on.Nav19 * 0.0002,'vshift',vshift);
            strct.adendOMLout.Nav19 = struct('gbar',  on.Nav19 * 0.0002,'vshift',vshift);
        end
        %% A-type Potassium
        
        % Kv1.1
        if on.Kv11 > 0
            kfac = on.Kv11;
            vshift = 0;
            gbdef = 0.002;
            strct.axon.Kv11 = struct('gkbar',gbdef*kfac,'vshift',vshift);%,'mk',-8.1,'md',-32);
            strct.axonh.Kv11 = struct('gkbar',gbdef*kfac,'vshift',vshift);%,'mk',-8.1,'md',-32);
        end
        
        % Kv1.4
        if on.Kv14 > 0 % Protein nearly only in axon!
            vshift = 0;
            strct.axon.Kv14 = struct('gkbar',0.008*on.Kv14,'vshift',vshift); %  scale_i von 3 möglich um inactivatin schneller zu machen wegen kvbeta1.1
            strct.axonh.Kv14 = struct('gkbar',0.008*on.Kv14,'vshift',vshift);
        end
        
        %Kv2.1
        if on.Kv21 > 0
            strct.soma.Kv21 = struct('gkbar',on.Kv21 * 0.01);
        end
        
        % Kv3.3 und Kv3.4 Protein nearly only in axon!
        if on.Kv33 > 0
            strct.axon.Kv33 = struct('gkbar',on.Kv33 * 0.002  );
            strct.axonh.Kv33 = struct('gkbar',on.Kv33 * 0.02 );
        end
        if on.Kv34 > 0
            if ~isempty(strfind(options,'-k'))
                scale_a = 2.5;
                ksl= 0.75;
            else
                scale_a = 4;%2.5;
                ksl= 0.5;%0.75;
            end
            strct.axon.Kv34 = struct('gkbar',on.Kv34 * 0.0049 ,'scale_a',scale_a,'ksl',ksl);  % 0.00325
            strct.axonh.Kv34 = struct('gkbar',on.Kv34 * 0.0197 ,'scale_a',scale_a,'ksl',ksl); % 0.0325
        end
        
        %Kv4.2
        if on.Kv42 > 0
            kafac = on.Kv42 * 1;
            gdist = [0.0005,0.001,0.001,0.001]; % GCL IML MML OML
            %                 gdist = [0.0002,0.001,0.003,0.003]; % GCL IML MML OML
            switch Kv4model{Kv4usethismodel}
                case 'Amarillo'
                    switch Kv4composition
                        case 1
                            % Kv4.2 alone
                            strct.GCL.Kv42b = struct('gkbar',kafac * gdist(1),'a0', 1.589,'za',0.64,'b0', 0.0184,'zb', -1.31,'c0',6.668,'zc', 0.15,'d0', 2.381,'zd', -1.21,'e0',0.503,'ze',0.07,'f0',0.174,'zf',-0.25,'kci', 0.047,'kic', 0.00003,'kappa1', 0.229,'lambda1', 0.151,'kappa2',0.0487,'lambda2',0.0065);
                            strct.adendIML.Kv42b = struct('gkbar',kafac * gdist(2),'a0', 1.589,'za',0.64,'b0', 0.0184,'zb', -1.31,'c0',6.668,'zc', 0.15,'d0', 2.381,'zd', -1.21,'e0',0.503,'ze',0.07,'f0',0.174,'zf',-0.25,'kci', 0.047,'kic', 0.00003,'kappa1', 0.229,'lambda1', 0.151,'kappa2',0.0487,'lambda2',0.0065);
                            strct.adendMML.Kv42b = struct('gkbar',kafac * gdist(3),'a0', 1.589,'za',0.64,'b0', 0.0184,'zb', -1.31,'c0',6.668,'zc', 0.15,'d0', 2.381,'zd', -1.21,'e0',0.503,'ze',0.07,'f0',0.174,'zf',-0.25,'kci', 0.047,'kic', 0.00003,'kappa1', 0.229,'lambda1', 0.151,'kappa2',0.0487,'lambda2',0.0065);
                            strct.adendOML.Kv42b = struct('gkbar',kafac * gdist(4),'a0', 1.589,'za',0.64,'b0', 0.0184,'zb', -1.31,'c0',6.668,'zc', 0.15,'d0', 2.381,'zd', -1.21,'e0',0.503,'ze',0.07,'f0',0.174,'zf',-0.25,'kci', 0.047,'kic', 0.00003,'kappa1', 0.229,'lambda1', 0.151,'kappa2',0.0487,'lambda2',0.0065);
                        case 3
                            % + DPPX
                            strct.GCL.Kv42b = struct('gkbar',kafac * gdist(1),'a0',3.282,'za',0.64,'b0', 0.01,'zb', -1.31,'c0',7.22,'zc', 0.15,'d0', 0.588,'zd', -1.21,'e0',1.101,'ze',0.07,'f0',0.673,'zf',-0.25,'kci', 0.072,'kic', 0.00046,'kappa1', 0.108,'lambda1', 0.034,'kappa2',0.0206,'lambda2',0.0086);
                            strct.adendIML.Kv42b = struct('gkbar',kafac * gdist(2),'a0',3.282,'za',0.64,'b0', 0.01,'zb', -1.31,'c0',7.22,'zc', 0.15,'d0', 0.588,'zd', -1.21,'e0',1.101,'ze',0.07,'f0',0.673,'zf',-0.25,'kci', 0.072,'kic', 0.00046,'kappa1', 0.108,'lambda1', 0.034,'kappa2',0.0206,'lambda2',0.0086);
                            strct.adendMML.Kv42b = struct('gkbar',kafac * gdist(3),'a0',3.282,'za',0.64,'b0', 0.01,'zb', -1.31,'c0',7.22,'zc', 0.15,'d0', 0.588,'zd', -1.21,'e0',1.101,'ze',0.07,'f0',0.673,'zf',-0.25,'kci', 0.072,'kic', 0.00046,'kappa1', 0.108,'lambda1', 0.034,'kappa2',0.0206,'lambda2',0.0086);
                            strct.adendOML.Kv42b = struct('gkbar',kafac * gdist(4),'a0',3.282,'za',0.64,'b0', 0.01,'zb', -1.31,'c0',7.22,'zc', 0.15,'d0', 0.588,'zd', -1.21,'e0',1.101,'ze',0.07,'f0',0.673,'zf',-0.25,'kci', 0.072,'kic', 0.00046,'kappa1', 0.108,'lambda1', 0.034,'kappa2',0.0206,'lambda2',0.0086);
                        case 2
                            % + KChIP1
                            strct.GCL.Kv42b = struct('gkbar',kafac * gdist(1),'a0',1.443,'za',0.64,'b0', 0.0019,'zb', -1.31,'c0',5.556,'zc', 0.15,'d0', 0.155,'zd', -1.21,'e0',0.559,'ze',0.07,'f0',0.621,'zf',-0.25,'kci', 0.029,'kic', 0.00014,'kappa1', 0,'lambda1', 100,'kappa2',0,'lambda2',100);
                            strct.adendIML.Kv42b = struct('gkbar',kafac * gdist(2),'a0',1.443,'za',0.64,'b0', 0.0019,'zb', -1.31,'c0',5.556,'zc', 0.15,'d0', 0.155,'zd', -1.21,'e0',0.559,'ze',0.07,'f0',0.621,'zf',-0.25,'kci', 0.029,'kic', 0.00014,'kappa1', 0,'lambda1', 100,'kappa2',0,'lambda2',100);
                            strct.adendMML.Kv42b = struct('gkbar',kafac * gdist(3),'a0',1.443,'za',0.64,'b0', 0.0019,'zb', -1.31,'c0',5.556,'zc', 0.15,'d0', 0.155,'zd', -1.21,'e0',0.559,'ze',0.07,'f0',0.621,'zf',-0.25,'kci', 0.029,'kic', 0.00014,'kappa1', 0,'lambda1', 100,'kappa2',0,'lambda2',100);
                            strct.adendOML.Kv42b = struct('gkbar',kafac * gdist(4),'a0',1.443,'za',0.64,'b0', 0.0019,'zb', -1.31,'c0',5.556,'zc', 0.15,'d0', 0.155,'zd', -1.21,'e0',0.559,'ze',0.07,'f0',0.621,'zf',-0.25,'kci', 0.029,'kic', 0.00014,'kappa1', 0,'lambda1', 100,'kappa2',0,'lambda2',100);
                            
                        case 4
                            % + KChIP1 and DPPX (ternary)
                            strct.GCL.Kv42b = struct('gkbar',kafac * gdist(1),'a0',2.577,'za',0.64,'b0', 0.0028,'zb', -1.31,'c0',4.318,'zc', 0.15,'d0', 0.38,'zd', -1.21,'e0',0.466,'ze',0.07,'f0',0.277,'zf',-0.25,'kci', 0.054,'kic', 0.00044,'kappa1', 0,'lambda1', 100,'kappa2',0,'lambda2',100);
                            strct.adendIML.Kv42b = struct('gkbar',kafac * gdist(2),'a0',2.577,'za',0.64,'b0', 0.0028,'zb', -1.31,'c0',4.318,'zc', 0.15,'d0', 0.38,'zd', -1.21,'e0',0.466,'ze',0.07,'f0',0.277,'zf',-0.25,'kci', 0.054,'kic', 0.00044,'kappa1', 0,'lambda1', 100,'kappa2',0,'lambda2',100);
                            strct.adendMML.Kv42b = struct('gkbar',kafac * gdist(3),'a0',2.577,'za',0.64,'b0', 0.0028,'zb', -1.31,'c0',4.318,'zc', 0.15,'d0', 0.38,'zd', -1.21,'e0',0.466,'ze',0.07,'f0',0.277,'zf',-0.25,'kci', 0.054,'kic', 0.00044,'kappa1', 0,'lambda1', 100,'kappa2',0,'lambda2',100);
                            strct.adendOML.Kv42b = struct('gkbar',kafac * gdist(4),'a0',2.577,'za',0.64,'b0', 0.0028,'zb', -1.31,'c0',4.318,'zc', 0.15,'d0', 0.38,'zd', -1.21,'e0',0.466,'ze',0.07,'f0',0.277,'zf',-0.25,'kci', 0.054,'kic', 0.00044,'kappa1', 0,'lambda1', 100,'kappa2',0,'lambda2',100);
                    end
                case 'Barghaan'
                    vshift = -15-5; % don't ask why but activation is ~ 15 mV right shifted in the model compared to the measurements (Suppl Fig1 Barghaan 2008), additional -5mV shift due to uncorrected LJP of 4.9
                    kafac = kafac * 10; %optional raise this..
                    switch Kv4composition
                        case 1
                            strct.GCL.Kv42 = struct('gkbar',kafac * gdist(1),'a0', 0.175,'za', 2.7,'b0', 0.003598,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.267,'zoc', -0.047,'kci', 0.02392,'kic', 0.000037,'koi', 0.194,'kio', 0.03686,'vshift',vshift);
                            strct.SGCL.Kv42 = struct('gkbar',kafac * gdist(1),'a0', 0.175,'za', 2.7,'b0', 0.003598,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.267,'zoc', -0.047,'kci', 0.02392,'kic', 0.000037,'koi', 0.194,'kio', 0.03686,'vshift',vshift);
                            strct.adendIML.Kv42 = struct('gkbar',kafac * gdist(2),'a0', 0.175,'za', 2.7,'b0', 0.003598,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.267,'zoc', -0.047,'kci', 0.02392,'kic', 0.000037,'koi', 0.194,'kio', 0.03686,'vshift',vshift);
                            strct.adendMML.Kv42 = struct('gkbar',kafac * gdist(3),'a0', 0.175,'za', 2.7,'b0', 0.003598,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.267,'zoc', -0.047,'kci', 0.02392,'kic', 0.000037,'koi', 0.194,'kio', 0.03686,'vshift',vshift);
                            strct.adendOML.Kv42 = struct('gkbar',kafac * gdist(4),'a0', 0.175,'za', 2.7,'b0', 0.003598,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.267,'zoc', -0.047,'kci', 0.02392,'kic', 0.000037,'koi', 0.194,'kio', 0.03686,'vshift',vshift);
                        case 2
                            strct.GCL.Kv42 = struct('gkbar',kafac * gdist(1),'a0', 0.175,'za', 2.7,'b0', 0.01947,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.670,'zoc', -0.047,'kci', 0.03061,'kic', 0.00018,'koi', 0.04308,'kio', 0.1099,'vshift',vshift);
                            strct.SGCL.Kv42 = struct('gkbar',kafac * gdist(1),'a0', 0.175,'za', 2.7,'b0', 0.01947,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.670,'zoc', -0.047,'kci', 0.03061,'kic', 0.00018,'koi', 0.04308,'kio', 0.1099,'vshift',vshift);
                            strct.adendIML.Kv42 = struct('gkbar',kafac * gdist(2),'a0', 0.175,'za', 2.7,'b0', 0.01947,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.670,'zoc', -0.047,'kci', 0.03061,'kic', 0.00018,'koi', 0.04308,'kio', 0.1099,'vshift',vshift);
                            strct.adendMML.Kv42 = struct('gkbar',kafac * gdist(3),'a0', 0.175,'za', 2.7,'b0', 0.01947,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.670,'zoc', -0.047,'kci', 0.03061,'kic', 0.00018,'koi', 0.04308,'kio', 0.1099,'vshift',vshift);
                            strct.adendOML.Kv42 = struct('gkbar',kafac * gdist(4),'a0', 0.175,'za', 2.7,'b0', 0.01947,'zb', -1.742,'kco0', 0.347,'zco', 0.185,'koc0', 1.670,'zoc', -0.047,'kci', 0.03061,'kic', 0.00018,'koi', 0.04308,'kio', 0.1099,'vshift',vshift);
                        case 3
                            strct.GCL.Kv42 = struct('gkbar',kafac * gdist(1),'a0', 0.416,'za', 1.1,'b0', 0.009,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.267,'zoc', 0,'kci', 0.03807,'kic', 0.00011,'koi', 0.3,'kio', 0.01424,'vshift',vshift);
                            strct.SGCL.Kv42 = struct('gkbar',kafac * gdist(1),'a0', 0.416,'za', 1.1,'b0', 0.009,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.267,'zoc', 0,'kci', 0.03807,'kic', 0.00011,'koi', 0.3,'kio', 0.01424,'vshift',vshift);
                            strct.adendIML.Kv42 = struct('gkbar',kafac * gdist(2),'a0', 0.416,'za', 1.1,'b0', 0.009,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.267,'zoc', 0,'kci', 0.03807,'kic', 0.00011,'koi', 0.3,'kio', 0.01424,'vshift',vshift);
                            strct.adendMML.Kv42 = struct('gkbar',kafac * gdist(3),'a0', 0.416,'za', 1.1,'b0', 0.009,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.267,'zoc', 0,'kci', 0.03807,'kic', 0.00011,'koi', 0.3,'kio', 0.01424,'vshift',vshift);
                            strct.adendOML.Kv42 = struct('gkbar',kafac * gdist(4),'a0', 0.416,'za', 1.1,'b0', 0.009,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.267,'zoc', 0,'kci', 0.03807,'kic', 0.00011,'koi', 0.3,'kio', 0.01424,'vshift',vshift);
                            
                        case 4
                            strct.GCL.Kv42 = struct('gkbar',kafac * gdist(1),'a0', 0.416,'za', 1.1,'b0', 0.0486,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.670,'zoc', 0,'kci', 0.04873,'kic', 0.000537,'koi', 0.0669,'kio', 0.04246,'vshift',vshift);
                            strct.adendIML.Kv42 = struct('gkbar',kafac * gdist(2),'a0', 0.416,'za', 1.1,'b0', 0.0486,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.670,'zoc', 0,'kci', 0.04873,'kic', 0.000537,'koi', 0.0669,'kio', 0.04246,'vshift',vshift);
                            strct.adendMML.Kv42 = struct('gkbar',kafac * gdist(3),'a0', 0.416,'za', 1.1,'b0', 0.0486,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.670,'zoc', 0,'kci', 0.04873,'kic', 0.000537,'koi', 0.0669,'kio', 0.04246,'vshift',vshift);
                            strct.adendOML.Kv42 = struct('gkbar',kafac * gdist(4),'a0', 0.416,'za', 1.1,'b0', 0.0486,'zb', -1.556,'kco0', 0.347,'zco', 0,'koc0', 1.670,'zoc', 0,'kci', 0.04873,'kic', 0.000537,'koi', 0.0669,'kio', 0.04246,'vshift',vshift);
                    end
            end
        end
        
        %Kv7.2/3 (M current)
        if on.Kv72u3 > 0
            Dtaumult = 6;
            tau0mult = 0.2;
            %                 on.Kv72u3 = on.Kv72u3 * 2
            strct.axonh.Kv723 = struct('gkbar',on.Kv72u3 *   0.0005 *2,'Dtaumult1',Dtaumult,'Dtaumult2',Dtaumult,'tau0mult',tau0mult);
            strct.axon.Kv723 = struct('gkbar',on.Kv72u3 *   0.0001 *2,'Dtaumult1',Dtaumult,'Dtaumult2',Dtaumult,'tau0mult',tau0mult);
        end
        
        %AH Calcium channels
        if on.TCa > 0 || on.NCa > 0 || on.LCa > 0

            strct.all.Cabuffer = struct('brat',50,'tau',Ca_tau ); 
            strct.axon.Cabuffer = struct('brat',10,'tau',43 ); 
            strct.soma.Cabuffer = struct('brat',200,'tau',Ca_tau );
            strct.all.Cabuffer.depth = 0.05; %µm   
            
            tfac = on.TCa ;
            fac = 10 * 1.5 ;  
            fac3 = fac * 5 / 6;
            on.LCa = on.LCa /30; 
            on.NCa = on.NCa *3;
            kf = 0.0005; % = 0.5µM = 500 nM
            vdi12 = 1;  % voltage dependent inactivation of the channel (0-1)
            vdi13 = 0.85;   % voltage dependent inactivation of the channel (0-1)
            strct.axon.Cav22 = struct('gbar', on.NCa * 0.0001 * 5 / fac,'hTau',80);
            strct.axonh.Cav22 = struct('gbar', on.NCa * 0.0001 * 5/ fac,'hTau',80);
            strct.soma.Cav22 = struct('gbar', on.NCa * 0.0015 *2 / fac,'hTau',80);
            strct.SGCL.Cav22 = struct('gbar', on.NCa * 0.0001 *5 / fac,'hTau',80);
            strct.GCL.Cav22 = struct('gbar', on.NCa * 0.0001 *5 / fac);
            strct.adendIML.Cav22 = struct('gbar', on.NCa * 0.0001 *5 / fac,'hTau',80);
            strct.adendMML.Cav22 = struct('gbar', on.NCa * 0.0001 *5 / fac,'hTau',80);
            strct.adendOML.Cav22 = struct('gbar', on.NCa * 0.0001 *5 / fac,'hTau',80);
            
            strct.axonh.Cav12 = struct('gbar', on.LCa * 0.0005 * 3/2 / fac,'kf',kf,'VDI',vdi12);
            strct.soma.Cav12 = struct('gbar', on.LCa  * 0.0015 / fac,'kf',kf,'VDI',vdi12);
            strct.SGCL.Cav12 = struct('gbar', on.LCa  * 0.0005 * 3/2 / fac,'kf',kf,'VDI',vdi12);
            strct.GCL.Cav12 = struct('gbar', on.LCa  * 0.0005 * 3/2 / fac,'kf',kf,'VDI',vdi12);
            strct.adendIML.Cav12 = struct('gbar', on.LCa  *  0.0005 * 3/2 * 4 / fac,'kf',kf,'VDI',vdi12);
            strct.adendMML.Cav12 = struct('gbar', on.LCa  * 0.0005 * 3/2 * 4 / fac,'kf',kf,'VDI',vdi12);
            strct.adendOML.Cav12 = struct('gbar', on.LCa  *  0.0005 * 3/2 * 4 / fac,'kf',kf,'VDI',vdi12);
            
            strct.axon.Cav13 = struct('gbar', on.LCa * 0.0005 / 2 / fac3,'kf',kf,'VDI',vdi13);
            strct.axonh.Cav13 = struct('gbar', on.LCa * 0.0005 / fac3,'kf',kf,'VDI',vdi13);
            strct.soma.Cav13 = struct('gbar', on.LCa  * 0.0005 * 2 / fac3,'kf',kf,'VDI',vdi13);
            strct.SGCL.Cav13 = struct('gbar', on.LCa  * 0.0005 /2  / fac3,'kf',kf,'VDI',vdi13);
            strct.GCL.Cav13 = struct('gbar', on.LCa  * 0.0005 /2 / fac3,'kf',kf,'VDI',vdi13);
            strct.adendIML.Cav13 = struct('gbar', on.LCa  *  0.0005 * 4/2 / 2 / fac3,'kf',kf,'VDI',vdi13);
            strct.adendMML.Cav13 = struct('gbar', on.LCa  * 0.0005 * 4/2 / 2 / fac3,'kf',kf,'VDI',vdi13);
            strct.adendOML.Cav13 = struct('gbar', on.LCa  *  0.0005 * 4/2 / 2 / fac3,'kf',kf,'VDI',vdi13);
            %
            strct.axon.Cav32 = struct('gbar', 0.0005/1.375 *tfac / fac);% no T-type in axon bouton (Li,Bischofberger 2007) but Martinello finds some immunogold...
            strct.axonh.Cav32 = struct('gbar', 0.0005/1.375  *tfac / fac);
            strct.soma.Cav32 = struct('gbar', 0.001 *tfac / fac);
            strct.SGCL.Cav32 = struct('gbar', 0.0002 * 2.5 *2 *tfac  / fac);
            strct.GCL.Cav32 = struct('gbar', 0.0002 * 2.5 *2 *tfac / fac);
            strct.adendIML.Cav32 = struct('gbar', 0.0002 * 2.5 *2 *tfac / fac);
            strct.adendMML.Cav32 = struct('gbar', 0.0002 * 2.5 *2 *tfac / fac);
            strct.adendOML.Cav32 = struct('gbar', 0.0002 * 2.5 *2 *tfac / fac);
        end
        
        if on.BK > 0
            gabk =  13 * 0.0003 * on.BK ;

            if ~isempty(strfind(options,'-n'))
                base = 4;
                diffsoma = 1.5;%0.3;
                gak = gabk * 4;
                
                gabk_soma = gabk / 4;
                gak_soma = gabk * 4 / 4;
            else
                base = 2;
                diffsoma = 1;%0.3;
                gabk_soma = gabk / 2;
                gak = 0;%gak / 2 ;
                gak_soma = 0;
            end

            % BK: 0.05 - 0.3 S/cm²! (from model)
            strct.axon.BK = struct('gabkbar',  gabk,'gakbar', gak,'diff',100000,'base',base);  % no Calcium channel - BK clusters in axon, this is why high diff value
            strct.axonh.BK = struct('gabkbar', gabk,'gakbar', gak,'diff',100000,'base',base); % axonh i dont know, but made it same as axon...
            % ab is slow bk type I, a is fast bk type II
            strct.soma.BK = struct('gabkbar', gabk_soma,'gakbar', gak_soma,'diff',diffsoma,'base',base);  % also more gak because there might be heteromeric BKs with less than 4 ß4 subunits
        end
        if on.SK > 0
            Q10 = 5;
            skfac = on.SK / 24; 
            dif = 3; 
            strct.soma.SK2 = struct('gkbar',0.004 / 20 * skfac,'diff',dif,'Q10',Q10,'fac',2.5,'diro2',0.1,'invc3',0.09,'invc1',0.16*2,'invc2',0.16*2,'dirc4',320); % this is not really SK2 but SK1
            strct.axonh.SK2 = struct('gkbar',0.02* skfac,'diff',dif,'Q10',Q10,'fac',2.5,'diro2',0.1,'invc3',0.09,'invc1',0.16*2,'invc2',0.16*2,'dirc4',320); % this is not really SK2 but SK1
            strct.axon.SK2 = struct('gkbar',0.003 * skfac,'diff',dif,'Q10',Q10,'fac',2.5,'diro2',0.1,'invc3',0.09,'invc1',0.16*2,'invc2',0.16*2,'dirc4',320); % this is not really SK2 but SK1/3
            strct.SGCL.SK2 = struct('gkbar',0.0002 * skfac,'diff',dif,'Q10',Q10,'fac',2.5,'diro2',0.1,'invc3',0.09,'invc1',0.16*2,'invc2',0.16*2,'dirc4',320);
            strct.GCL.SK2 = struct('gkbar',0.0004 * skfac,'diff',dif,'Q10',Q10,'fac',2.5,'diro2',0.1,'invc3',0.09,'invc1',0.16*2,'invc2',0.16*2,'dirc4',320);
            
            strct.adendIML.SK2 = struct('gkbar',0.00105 * skfac,'diff',dif,'Q10',Q10,'fac',2.5,'diro2',0.1,'invc3',0.09,'invc1',0.08*4,'invc2',0.08*4,'dirc4',80*4);  % ~ 35 mal weniger SK3 in soma als in IML/OML (Ballesteros-Merino 2014)
            strct.adendMML.SK2 = struct('gkbar',0.00105 * skfac,'diff',dif,'Q10',Q10,'fac',2.5,'diro2',0.1,'invc3',0.09,'invc1',0.16*2,'invc2',0.16*2,'dirc4',320);
            strct.adendOML.SK2 = struct('gkbar',0.00105 * skfac,'diff',dif,'Q10',Q10,'fac',2.5,'diro2',0.1,'invc3',0.09,'invc1',0.16*2,'invc2',0.16*2,'dirc4',320);
        end

        strct.all.k_ion = struct('ek',-103);      % calculated from concentrations in Mongiat paper (SH07 very similar)
        strct.all.na_ion = struct('ena',87);     % calculated from concentrations in Mongiat paper (SH07 very similar)
        strct.all.ca_ion = struct('cao0',2,'cai0',0.00007);   % calculated from concentrations in Mongiat paper (SH07 very similar)
        
    else
        %already with corrected Ca buffer model...
        strct.axon.ichan3 = struct('gnabar',0.21,'gkfbar',0.028,'gksbar',0,'gkabar',0.004);
        strct.axonh.ichan3 = struct('gnabar',0.21,'gkfbar',0.028,'gksbar',0,'gkabar',0.004);
        strct.soma.ichan3 = struct('gnabar',0.12,'gkfbar',0.016,'gksbar',0.003,'gkabar',0.012);
        strct.soma.Caold = struct('gtcabar',0.00015,'gncabar',0.002,'glcabar',0.01);
        strct.soma.CadepK = struct('gbkbar',0.0003,'gskbar',0.0005);
        strct.GCL.ichan3 = struct('gnabar',0.018,'gkfbar',0.004,'gksbar',0.003,'gkabar',0);
        strct.SGCL.ichan3 = struct('gnabar',0.018,'gkfbar',0.004,'gksbar',0.003,'gkabar',0);
        strct.GCL.Caold = struct('gtcabar',0.0003,'gncabar',0.003,'glcabar',0.015);
        strct.GCL.CadepK = struct('gbkbar',0.0003,'gskbar',0.0002);
        strct.SGCL.Caold = struct('gtcabar',0.0003,'gncabar',0.003,'glcabar',0.015);
        strct.SGCL.CadepK = struct('gbkbar',0.0003,'gskbar',0.0002);
        strct.adendIML.ichan3 = struct('gnabar',0.013,'gkfbar',0.004,'gksbar',0.003,'gkabar',0);
        strct.adendIML.Caold = struct('gtcabar',0.001,'gncabar',0.001,'glcabar',0.015);
        strct.adendIML.CadepK = struct('gbkbar',0.0005,'gskbar',0.0001);
        strct.adendMML.ichan3 = struct('gnabar',0.008,'gkfbar',0.001,'gksbar',0.003,'gkabar',0);
        strct.adendMML.Caold = struct('gtcabar',0.002,'gncabar',0.001,'glcabar',0.001);
        strct.adendMML.CadepK = struct('gbkbar',0.0012,'gskbar',0);
        strct.adendOML.ichan3 = struct('gnabar',0,'gkfbar',0.001,'gksbar',0.004,'gkabar',0);
        strct.adendOML.Caold = struct('gtcabar',0.002,'gncabar',0.001,'glcabar',0);
        strct.adendOML.CadepK = struct('gbkbar',0.0012,'gskbar',0);  
        strct.all.k_ion = struct('ek',-85);     % calculated from concentrations in Mongiat paper (SH07 very similar)
        strct.all.na_ion = struct('ena',45);   %  % calculated from concentrations in Mongiat paper (SH07 very similar)
        strct.all.ca_ion = struct('cao0',2,'cai0',0.00007);    % calculated from concentrations in Mongiat paper (SH07 very similar)
    end
end