function aGC_synstim(neuron,tree,type,ostruct,targetfolder_data)
if nargin < 5
    targetfolder_data = 'D:/EphysModel';
end

if isstruct(neuron)
    neuron = {neuron};
else
    warndlg('why is neuron already a cell structure?')
end

% initial syn parameters...are changed later
%
switch ostruct.synmode
    case 1
        ppweight = ones(1,numel(tree)) * 0.00065; % SH07 = 0.1-1nS ....%weight =  µS
        nsyn = 30;
    case 2
        if ostruct.newborn
            ppweight = ones(1,numel(tree)) * 0.00065; % SH07 = 0.1-1nS ....%weight =  µS
            nsyn = 15;
        else
            ppweight = ones(1,numel(tree)) * 0.00065; % SH07 = 0.1-1nS ....%weight =  µS
            nsyn = 30;
        end
    case 3
        nsyn = 30;
        ppweight = 0.0003;% SH07 = 0.1-1nS ....%weight =  µS
end

freq = NaN;
dd0 = NaN;
dt0 = NaN;

neuron{1}.params.dt = 0.05;
neuron{1}.params.cvode = 0;

recnode = cell(numel(tree),1);
if ~isempty(strfind(type,'spatial'))
    dd0 = 0:10:100;  % µm !!
    thesesynidTags = cell(numel(dd0),numel(tree));
else
    thesesynidTags = recnode;
end
s = RandStream.create('mt19937ar','Seed',1204); % be sure to always hit the same "random" nodes in order to assure comparability

for t = 1:numel(tree)
    plen = Pvec_tree(tree{t});
    recnode{t} = [1,find(abs(plen-100)<0.5 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon')),1),find(abs(plen-200)<0.5 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon')),1)];
    neuron{1}.record{t}.cell = struct('node',recnode{t},'record',{'cai','v'});
    
    switch type
        
        case 'spatial'
            dirs = {'X','Y','Z'};
            % evt sollte ich anfang MML nehmen und von da ab hoch gehen...
%             ppweight = 0.0005;% SH07 = 0.1-1nS ....%weight =  µS
            idpar = idpar_tree(tree{t});
            synids1 = find(tree{t}.R == find(strcmp(tree{t}.rnames,'adendMML')) & tree{t}.R(idpar) == find(strcmp(tree{t}.rnames,'adendIML'))); %get all branch starts in MML
            for n = 1:numel(dd0)
                [~, ind] = min([std(tree{t}.X(synids1)),std(tree{t}.Y(synids1)),std(tree{t}.Z(synids1))]); % find direction of tree layering
                synids2 = find(abs(tree{t}.(dirs{ind}) - mean(tree{t}.(dirs{ind})(synids1)) - dd0(n)) < 1 & (tree{t}.R == find(strcmp(tree{t}.rnames,'adendMML')) | tree{t}.R == find(strcmp(tree{t}.rnames,'adendOML')))); % find nodes at a distance dd0 from the IML/MML border
                synids2 = setdiff(synids2,idpar(synids2));  % delete all direct parent nodes (due to rough distance search)
                thesesynids = [synids1(1:min(ceil(nsyn/2),numel(synids1)));synids1(randi(s,numel(synids1),ceil(nsyn/2)-numel(synids1),1));synids2(1:min(floor(nsyn/2),numel(synids2)));synids2(randi(s,numel(synids2),floor(nsyn/2)-numel(synids2),1))]';  % ensure that each branch is there at least once
                [neuron{n}.pp{t}.Exp2Syn , thesesynidTags{n,t}] = addExp2Syn(t,thesesynids);
            end
        case 'spatial2'
                        dirs = {'X','Y','Z'};
            % evt sollte ich anfang MML nehmen und von da ab hoch gehen...
%             ppweight = 0.0005;% SH07 = 0.1-1nS ....%weight =  µS
            idpar = idpar_tree(tree{t});

            synids1 = find(tree{t}.R == find(strcmp(tree{t}.rnames,'adendMML')) & tree{t}.R(idpar) == find(strcmp(tree{t}.rnames,'adendIML'))); %get all branch starts in MML
            ipar = ipar_tree(tree{t});
            [~, ind] = min([std(tree{t}.X(synids1)),std(tree{t}.Y(synids1)),std(tree{t}.Z(synids1))]); % find direction of tree layering
            [~,thisone] = max(tree{t}.(dirs{ind}));
            [~,~,thisone] = intersect(ipar(thisone,:),synids1);
            for n = 1:numel(dd0)
                synids2 = find(any(ipar == synids1(thisone),2) & abs(tree{t}.(dirs{ind}) - mean(tree{t}.(dirs{ind})(synids1)) - dd0(n)) < 1);
                synids2 = setdiff(synids2,idpar(synids2));  % delete all direct parent nodes (due to rough distance search)
                thesesynids = [synids1(thisone);synids2(1)]';  % ensure that each branch is there at least once
                [neuron{n}.pp{t}.Exp2Syn , thesesynidTags{n,t}] = addExp2Syn(t,thesesynids);
            end
        case 'test'
            synids1 = find(abs(plen-150)<1 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon')));
            thesesynids = synids1(randi(s,numel(synids1),nsyn,1));
            [neuron{1}.pp{t}.Exp2Syn , thesesynidTags{t}] = addExp2Syn(t,thesesynids);
        otherwise
            synids1 = find(tree{t}.R == find(strcmp(tree{t}.rnames,'adendMML')));  %
            synids2 = find(tree{t}.R == find(strcmp(tree{t}.rnames,'adendOML')));  %
            thesesynids = [synids1(randi(s,numel(synids1),ceil(nsyn/2),1));synids2(randi(s,numel(synids2),floor(nsyn/2),1))]; 
            [neuron{1}.pp{t}.Exp2Syn , thesesynidTags{t}] = addExp2Syn(t,thesesynids);
    end
end

if ostruct.synmode == 3
        if ~isempty(strfind(tree{1}.name,'SH07all2'))
                switch type
                    case 'temporal'
                        if ostruct.newborn
                            %                     ppweight = t2n_findSubthreshWeight(neuron,tree,ppweight,10);
                            ppweight = [0.3000    0.4500    0.3125    0.3000    0.3625    0.3000    0.3125    0.2125]*1e-3;
                            ppweight = ppweight * 0.9;
                        else
                            ppweight = [0.6500    1.0000    0.7500    0.6625    0.8500    0.5625    0.7125    0.5625]*1e-3;
                            ppweight = ppweight * 0.9;
                        end
                    case 'spatial'
                        if ostruct.newborn
                            %                     ppweight = t2n_findSubthreshWeight(neuron,tree,ppweight,10);
                            ppweight = [0.3000    0.4500    0.3125    0.3000    0.3125    0.3000    0.3000    0.2125]*1e-3;
                            ppweight = ppweight * 0.9;
                        else
                            ppweight = [0.6000    0.9625    0.7000    0.6000    0.7125    0.5500    0.6625    0.5000]*1e-3;
                            ppweight = ppweight * 0.9;
                        end
                end
        elseif ~isempty(strfind(tree{1}.name,'mouse_matGC_art'))
                dashgfdgd
        elseif ~isempty(strfind(tree{1}.name, 'Beining'))
                switch type
                    case 'temporal'
                        if ostruct.newborn
                            %ppweight = t2n_findSubthreshWeight(neuron,tree,ppweight,10);
                            ppweight = []*1e-3;
                        else
                            %ppweight = t2n_findSubthreshWeight(neuron,tree,ppweight,10);
                            ppweight = []*1e-3;
                        end
                    case 'spatial'
                        %ppweight = t2n_findSubthreshWeight(neuron,tree,ppweight,10);
                        ppweight = []*1e-3;
                end
        end
end

%% Test pulse protocol
switch type
    case 'test'
        indstim = numel(tree)+1;
        tree{indstim} = struct('artificial','NetStim','start',10,'interval',100,'number',1);
        if isfield(neuron{1},'con')
            neuron{1} = rmfield(neuron{1},'con');
        end
        for t = 1:numel(tree)-1
            neuron{1}.con(t) = struct('source',struct('cell',indstim,'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}),'weight',ppweight(t),'delay',0,'threshold',0.5);
        end
        
        neuron{1}.params.tstop = 50;
    case 'white'
        freq = 50; % Hz later -> MHz
        neuron{1}.params.cvode = 0;  % with white noise, no cvode recommended for high nsyn
        neuron{1}.params.tstop = 1000;
        indstim = numel(tree)+1:numel(tree)+nsyn;
        tree(indstim) = {struct('artificial','VecStim')};
        [spikeMat,tvec] = t2n_poissonSpikeGen(freq,neuron.params,numel(indstim));
        for in = 1:numel(indstim)
            
            
            neuron{1}.play{indstim(in)}.cell =  struct('node',1,'play','spike','times',tvec(spikeMat(in,:)));
            for t = 1:numel(tree)-numel(indstim)
                if ~isfield(neuron,'con')
                    neuron{1}.con(1) = struct('source',struct('cell',indstim(in),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}(in)),'weight',ppweight(t),'delay',0,'threshold',0.5);
                else
                    neuron{1}.con(end+1) = struct('source',struct('cell',indstim(in),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}(in)),'weight',ppweight(t),'delay',0,'threshold',0.5);
                end
            end
        end
    case 'regular'
        freq = [5,10,20,40,75,100]; % Hz later -> MHz

        neuron{1}.params.cvode = 0;  % with white noise, no cvode recommended for high nsyn
        neuron{1}.params.tstop = 1000;
        indstim = numel(tree)+1;
        tree(indstim) = {struct('artificial','VecStim')};
        for t = 1:numel(tree)-1
            if ~isfield(neuron{1},'con')
                neuron{1}.con(1) = struct('source',struct('cell',indstim,'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}),'weight',ppweight(t),'delay',0,'threshold',0.5);
            else
                neuron{1}.con(end+1) = struct('source',struct('cell',indstim,'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}),'weight',ppweight(t),'delay',0,'threshold',0.5);
            end
        end
        for f = 1:numel(freq)
            neuron{f}.play{indstim}.cell = struct('node',1,'play','spike','times',1000/freq(f):1000/freq(f):neuron{f}.params.tstop);
        end
        neuron = t2n_as(1,neuron);
    case 'temporal'
        freq = [10,20,40,75]; % Hz later -> MHz
            dt0 = -25:5:25; % ms
%             dt0 = -50:5:50;
        neuron{1}.params.cvode = 0;  % with white noise, no cvode recommended for high nsyn
        neuron{1}.params.tstop = 500;
        indstim = numel(tree)+1:numel(tree)+2;
        tree(indstim(1)) = {struct('artificial','VecStim')};
        tree(indstim(2)) = {struct('artificial','VecStim')};
        for t = 1:numel(tree)-2
            if ~isfield(neuron{1},'con')
                neuron{1}.con(1) = struct('source',struct('cell',indstim(1),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}(1:ceil(nsyn/2))),'weight',ppweight(t),'delay',0,'threshold',0.5);  % connect stim1 to MML syns
            else
                neuron{1}.con(end+1) = struct('source',struct('cell',indstim(1),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}(1:ceil(nsyn/2))),'weight',ppweight(t),'delay',0,'threshold',0.5); % connect stim1 to MML syns
            end
            neuron{1}.con(end+1) = struct('source',struct('cell',indstim(2),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}(ceil(nsyn/2)+1:end)),'weight',ppweight(t),'delay',0,'threshold',0.5); % connect stim2 to OML syns
        end
        for f = 1:numel(freq)
            for n = 1:numel(dt0)
                ind = (f-1)*numel(dt0) + n;
                neuron{ind}.play{indstim(1)}.cell = struct('node',1,'play','spike','times',(1000/freq(f):1000/freq(f):neuron{1}.params.tstop));
                neuron{ind}.play{indstim(2)}.cell = struct('node',1,'play','spike','times',(1000/freq(f):1000/freq(f):neuron{1}.params.tstop)+dt0(n));                    
            end
        end
        neuron = t2n_as(1,neuron);
    case 'spatial'

        freq = [10,20,40,75]; % Hz later -> MHz
        neuron{1}.params.cvode = 0;  % with white noise, no cvode recommended for high nsyn
        neuron{1}.params.tstop = 500;
        indstim = numel(tree)+1:numel(tree)+2;
        tree(indstim(1)) = {struct('artificial','VecStim')};
        for t = 1:numel(tree)-1
            for n = 1:numel(dd0)
                if ~isfield(neuron{n},'con') || ~isstruct(neuron{n}.con)
                    neuron{n}.con(1) = struct('source',struct('cell',indstim(1),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{n,t}),'weight',ppweight(t),'delay',0,'threshold',0.5);  % connect stim1 to MML syns %(1:ceil(nsyn/2))
                else
                    neuron{n}.con(end+1) = struct('source',struct('cell',indstim(1),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{n,t}),'weight',ppweight(t),'delay',0,'threshold',0.5); % connect stim1 to MML syns
                end
            end
        end
        for f = 1:numel(freq)
            for n = 1:numel(dd0)
                ind = (f-1)*numel(dd0) + n;
                neuron{ind}.play{indstim(1)}.cell = struct('node',1,'play','spike','times',(1000/freq(f):1000/freq(f):neuron{1}.params.tstop));
                if f > 1
                    neuron{ind} = t2n_as(n,neuron{ind});
                else
                    neuron{ind} = t2n_as(1,neuron{ind});
                end
            end
        end
    case 'spatial2'

        freq = 10;%,20,40,75]; % Hz later -> MHz
            ppweight = ones(1,numel(tree))*0.0065;
        neuron{1}.params.cvode = 0;  % with white noise, no cvode recommended for high nsyn
        neuron{1}.params.tstop = 500;
        indstim = numel(tree)+1:numel(tree)+2;
        tree(indstim(1)) = {struct('artificial','VecStim')};
        for t = 1:numel(tree)-1
            for n = 1:numel(dd0)
                if ~isfield(neuron{n},'con') || ~isstruct(neuron{n}.con)
                    neuron{n}.con(1) = struct('source',struct('cell',indstim(1),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{n,t}),'weight',ppweight(t),'delay',0,'threshold',0.5);  % connect stim1 to MML syns %(1:ceil(nsyn/2))
                else
                    neuron{n}.con(end+1) = struct('source',struct('cell',indstim(1),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{n,t}),'weight',ppweight(t),'delay',0,'threshold',0.5); % connect stim1 to MML syns
                end
            end
        end
        for f = 1:numel(freq)
            for n = 1:numel(dd0)
                ind = (f-1)*numel(dd0) + n;
                neuron{ind}.play{indstim(1)}.cell = struct('node',1,'play','spike','times',(1000/freq(f):1000/freq(f):neuron{1}.params.tstop));
                if f > 1
                    neuron{ind} = t2n_as(n,neuron{ind});
                else
                    neuron{ind} = t2n_as(1,neuron{ind});
                end
            end
        end
    case 'TBS' % TBS-protocol Ge 2007
        delay=200;
        interval1 = 10000;
        interval2 = 200;
        interval3 = 10;
        repeat1 = 4;
        repeat2 = 10;
        repeat3 = 10;
        
        block=delay:repeat3*interval3:repeat2*interval2+delay-1;
        times=[];
        for b=1:repeat1
            times = cat(2,times,block+interval1*(b-1));
        end
        
        indstim = numel(tree)+1:numel(tree)+4;
        
        tree{indstim(1)} = struct('artificial','NetStim','start',delay,'interval',interval1,'number',repeat1);
        tree{indstim(2)} = struct('artificial','NetStim','start',-1,'interval',interval2,'number',repeat2);
        tree{indstim(3)} = struct('artificial','NetStim','start',-1,'interval',interval3,'number',repeat3);
        tree{indstim(4)} = struct('artificial','NetStim','start',10,'interval',times(end)+300,'number',2);  % test stim at beginning and end of protocol
        
        
        if isfield(neuron{1},'con')
            neuron{1} =  rmfield(neuron{1},'con');
        end
        
        for t = 1:numel(tree)-numel(indstim)
            neuron{1}.con(2*(t-1)+1) = struct('source',struct('cell',indstim(3),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}),'weight',ppweight(t),'delay',0,'threshold',0.5);
            neuron{1}.con(2*t) = struct('source',struct('cell',indstim(4),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','tag',thesesynidTags{t}),'weight',ppweight(t),'delay',0,'threshold',0.5);
            neuron{1}.pp{t}.IClamp = struct('node',1,'del',-1e5,'dur',1e15,'times',[-100,times],'amp', [0,repmat([0.1,0],[1,repeat1*repeat2])]);  % paired 100pA current injection
        end
        
        neuron{1}.con(end+1) = struct('source',struct('cell',indstim(1),'watch','on'),'target',struct('cell',indstim(2)),'weight',1,'delay',0,'threshold',0.5);%,0.5,0,1});
        neuron{1}.con(end+1) = struct('source',struct('cell',indstim(2),'watch','on'),'target',struct('cell',indstim(3)),'weight',1,'delay',0,'threshold',0.5);%,0.5,0,1});

        neuron{1}.params.tstop = times(end)+400;
        neuron{1}.params.tstop = 2500;
end

[out, ~] = t2n(neuron,tree,'-w-q-d');
str = '';

if ostruct.newborn
    if ostruct.newborn == 2
        str = strcat(str,'_newborn_nopas');
    else
        str = strcat(str,'_newborn');
    end
end

if ~ (ostruct.vmodel >= 0)
    str = strcat(str,'_AH99');
end

switch ostruct.synmode 
    case 3
        str2 = '_subthresh';
    case 2
        str2 = '_newbornadj';
    otherwise
        str2 =  '_noadj';
end

switch type
    case 'temporal'
        str3 = dt0(end);
    case 'spatial'
        str3 = unique(diff(dd0));     
    otherwise
        str3='';
end

str = sprintf('%s_%gsyn_%gnS',str,nsyn,mean(ppweight)*1000);
save(fullfile(targetfolder_data,sprintf('Exp_%s_%s%s%s%d.mat',type,neuron{1}.experiment,str,str2,str3)),'out','neuron','tree','indstim','recnode','thesesynids','nsyn','ppweight','freq','dd0','dt0')
end

function [strct, tags] = addExp2Syn(t,theseSynIds)
    tags = string(compose('Tree%03d_syn%d',t,1:numel(theseSynIds)))';
    strct = struct('node',theseSynIds,'tag',tags,'tau1',0.2,'tau2',2.5,'e',0);%,'i',0.05,'e',0,'tau1',0.05,'tau2',2); tau values as SH07 %alt:tau1 1.5 tau2 5.5
end