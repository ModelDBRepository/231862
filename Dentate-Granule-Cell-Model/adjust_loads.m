function neuron = adjust_loads(neuron,tree,animal,ostruct)
% mode 1 adjust axon/soma conductances
% mode 2 adjust Ra

if nargin < 3
    animal = 'm';
end

region = {'axonh','soma','range'};
corresp_fields = {'Rho_AIS','Rho_soma'};
adjust_channels =  {'all'};{'na8st','Kv34','Kv14','BK'};%,'HKir'}; % only scaling for channels that contribute directly to spiking
% {'na8st','Kv34','Kv14','BK','SK2','Kv723','Cav12','Cav13','Cav22','Cav32'};

switch animal
    case 'm'
        Rho_exemplary = [100,19];%%%!!;  % Rho value for axon vs soma % 13 is too high as width goes up and 2nd amp down
    case 'r'
        if nargin < 3 || ostruct.ratadjust == 0
            Rho_exemplary = [100,19];  % Rho value for axon vs soma
        else
            Rho_exemplary = [100 19];[80 23];[68.4,27.2]; [40,6];  % Rho value for axon vs soma
        end
end
% adjust_channels = {'na8st','Kv33net','Kv34net','Kv14','Kv723','Ca2','Cav32','BK','SK2','Kv21'};%,'HKir'}; % no "passive" channel scaling
% noadjustchan = {'HKir','Kir21','pas'};% no "passive" channel scaling
for t = 1:numel(neuron.mech)
    
    adjust_region = region;
    
    if any(strcmp(adjust_region,'range'))  % get the index nodes for all regions that are adjusted...needed later for the ranging variables
        ind = false(numel(tree{t}.X),numel(adjust_region)-1);
        for f1 = 1:numel(adjust_region)-1
            ind(:,f1) = tree{t}.R == find(strcmp(tree{t}.rnames,adjust_region{f1})); % get index to all nodes of that region
        end
    end
    
    for f1 = 1:numel(adjust_region)
        if isfield(neuron.mech{t},adjust_region{f1})
            fields2 = fieldnames(neuron.mech{t}.(adjust_region{f1}));
            if ~any(strcmp(adjust_channels,'all'))
                fields2 = intersect(fields2,adjust_channels);%setdiff(fields2,noadjustchan);%intersect(fields2,adjust_channels);
            end
            for f2 = 1:numel(fields2)
                fields3 = fieldnames(neuron.mech{t}.(adjust_region{f1}).(fields2{f2}));  % get all parameter region to look for "bar"-words
                for f3 = 1:numel(fields3)
                    if ~isempty(strfind(fields3{f3},'bar'))  % look for "bar" word in parameter. if yes scale that parameter
                        if strcmp(adjust_region{f1},'range')
                            % more complex handling because the bar value
                            % is now for each node...
                            for ff1 = 1:numel(adjust_region)-1
                                neuron.mech{t}.(adjust_region{f1}).(fields2{f2}).(fields3{f3})(ind(:,ff1)) = neuron.mech{t}.(adjust_region{f1}).(fields2{f2}).(fields3{f3})(ind(:,ff1)) * tree{t}.(corresp_fields{ff1})/Rho_exemplary(ff1);  % scale gbar according to ratio this Rho vs the exemplary cell's Rho
                            end
                        else
%                             if strcmp(fields2{f2},'na8st')
                                neuron.mech{t}.(adjust_region{f1}).(fields2{f2}).(fields3{f3}) = neuron.mech{t}.(adjust_region{f1}).(fields2{f2}).(fields3{f3}) * tree{t}.(corresp_fields{f1})/Rho_exemplary(f1);  % scale gbar according to ratio this Rho vs the exemplary cell's Rho
%                             else
%                                 neuron.mech{t}.(adjust_region{f1}).(fields2{f2}).(fields3{f3}) = neuron.mech{t}.(adjust_region{f1}).(fields2{f2}).(fields3{f3}) * Rho_exemplary(f1) / tree{t}.(corresp_fields{f1});  % scale gbar according to ratio this Rho vs the exemplary cell's Rho
%                             end
                        end
                    end
                end
                
            end
        end
    end
    % dämliches ding für range....
    if isfield(neuron.mech{t},'range')
        fields2 = fieldnames(neuron.mech{t}.range);
        for f2 = 1:numel(fields2)
            fields3 = fieldnames(neuron.mech{t}.range.(fields2{f2}));
            for f3 = 1:numel(fields3)
                if ~isempty(strfind(fields3{f3},'bar'))
                    
                end
            end
        end
    end
end


if isfield(neuron,'experiment') && ~isempty(neuron.experiment)
    neuron.experiment = strcat(neuron.experiment,'_Loadadjusted');
else
    neuron.experiment = 'Loadadjusted';
end
