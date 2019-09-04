function [thisdata,steps,rate] = load_ephys(datanum,str,extractkir)

if nargin < 3
    extractkir = 0;
end

folder = 'raw data';  % folder where raw data was placed

switch datanum
    case 1
        dat = load(sprintf('%s%sMongiat_Mature_%s.mat',folder,filesep,str));
        thisdata = dat.data{1};
        rate = dat.rate;
%         steps = dat.csteps;
        if strcmp(str,'CClamp')
            thisdata = thisdata-5;  % don't know why but data is systematically shifted by +5 mV. You see it because it doesnt start at -70 mV
        end
%     case 2.21
%         dat = load(sprintf('%s%sMongiat_Young_21dpi_%s.mat',folder,filesep,str));
%         thisdata = dat.data{1};
%         rate = dat.rate;
% %         steps = dat.csteps;
%     case 2.25
%         dat = load(sprintf('%s%sMongiat_Young_25dpi_%s.mat',folder,filesep,str));
%         thisdata = dat.data{1};
%         rate = dat.rate;
% %         steps = dat.csteps;
    case 2.28
        dat = load(sprintf('%s%sMongiat_Young_28dpi_%s.mat',folder,filesep,str));
        thisdata = dat.data{1};
        rate = dat.rate;
%         steps = dat.csteps;
    case 2
        dat = load(sprintf('%s%sMongiat_Young_%s.mat',folder,filesep,str));
        thisdata = dat.data{1};
        rate = dat.rate;
%         steps = dat.csteps;
    case 3
        dat = load(sprintf('%s%sMongiat_BaCl_%s.mat',folder,filesep,str));
        if extractkir && strcmp(str,'VClamp')
            thisdata = dat.data{2}-dat.data{4};  % subtract the Ba insensitive current
        else
            thisdata = dat.data{2};
        end
        rate = dat.rate;
%         steps = dat.csteps;
    case 4
        dat = load(sprintf('%s%sMongiat_BaCl_%s.mat',folder,filesep,str));
        if extractkir
            thisdata = dat.data{1}-dat.data{3};  % subtract the Ba insensitive current
        else
            thisdata = dat.data{1};
        end
        rate = dat.rate;
%         steps = dat.csteps;
    case 5
        dat = load(sprintf('%s%sMongiat_BaCl_%s.mat',folder,filesep,str));
        thisdata = dat.data{4};
        rate = dat.rate;
%         steps = dat.csteps;
    case 6
        dat = load(sprintf('%s%sMongiat_BaCl_%s.mat',folder,filesep,str));
        thisdata = dat.data{3};
        rate = dat.rate;
%         steps = dat.csteps;
    case 7
        steps = (50:400)/1000;
        thisdata = NaN(1,1,8);

    case 0
        steps = [];
        thisdata = [];
        rate = [];
    otherwise
        error('datanumber not found')
end

switch str
    case 'VClamp'
        if datanum < 7
            steps = (-130:5:-40) - 12.1; %LJP corrected
        else
           steps = [];
        end
        
    case 'CClamp'
        if datanum < 7
            steps = (0:5:120)/1000;
            thisdata = thisdata - 12.1; %LJP corrected
        end
end
