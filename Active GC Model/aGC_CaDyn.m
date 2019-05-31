function aGC_CaDyn(neuron,tree,targetfolder_data,ostruct)

neuron.params.v_init = -85.4;

if ~isfield(ostruct,'cai')
    ostruct.cai = 'cai';
end
if ~isfield(ostruct,'cstep')
    ostruct.cstep = 1.3; % nA
end

neuron.params.tstop = 1000;
neuron.params.dt=0.025;
neuron.params.cvode = 1;
nodes = cell(numel(tree),1);
plen = nodes;
eucl = nodes;
plotcaivals = nodes;
CaNodes = nodes;
ipar = nodes;

hstep = t2n_findCurr(neuron,tree,neuron.params.v_init); %assuming a HP of xxx mV


for t = 1:numel(tree)
    
    plen{t} = Pvec_tree(tree{t});
    CaNodes{t} = {find(abs(plen{t}-150)<0.5 & tree{t}.R == find(strcmp(tree{t}.rnames,'axon'))),find(plen{t} == 0,1,'first'),find(abs(plen{t}-70)<0.5 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon'))),find(abs(plen{t}-150)<0.5 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon')))};
    if any(cellfun(@isempty,CaNodes{t}))
        fprintf('thresh distance of 0.5 for tree %d did not suffice. switched to thresh of 1\n',t)
        CaNodes{t} = {find(abs(plen{t}-150)<1 & tree{t}.R == find(strcmp(tree{t}.rnames,'axon'))),find(plen{t} == 0,1,'first'),find(abs(plen{t}-70)<1 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon'))),find(abs(plen{t}-150)<1 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon')))};
    end
    sprintf('Proximal tree %d: %s',t,strcat(tree{t}.rnames{unique(tree{t}.R(CaNodes{t}{3}))}))
    sprintf('Distal tree %d: %s',t,strcat(tree{t}.rnames{unique(tree{t}.R(CaNodes{t}{4}))}))
    
    ipar{t} = ipar_tree(tree{t});
    ipar{t} = ipar{t}(T_tree(tree{t}),:);  % only paths from termination points
    ipar{t}(ipar{t}==0) = 1;
    
    if ostruct.simple
        uipar = ipar(1,:);
        nodes{t} = uipar(1:find(uipar==1,1,'first'));
    else
        nodes{t} = unique(ipar{t});
    end
    
    % this part would have been done by t2n anyway, however to avoid
    % loading a lot of redundant values into Matlab, nodes are reduced to
    % the locations were NEURON actually calculates voltage here
    minterf = load(fullfile(pwd,'morphos','hocs',sprintf('%s_minterf.mat',tree{t}.NID)));
    minterf = t2n_makeNseg(tree{t},minterf.minterf,neuron.params,neuron.mech{t});
    inode = zeros(numel(nodes{t}),1);
    for in = 1:numel(nodes{t})
        inode(in) = find(minterf(:,1) == nodes{t}(in),1,'first');    %find the index of the node in minterf
    end
    [~,ia] = unique(minterf(inode,[2,4]),'rows');
    nodes{t} = sort(nodes{t}(ia));
    if ostruct.reduce  % reduce number of real recorded nodes to every third.
        nodes{t} = nodes{t}(1:3:end);
    end
    
    neuron.record{t}.cell(1:2) = struct('node',unique(cat(1,CaNodes{t}{:},nodes{t})),'record',{ostruct.cai,'v'});
    neuron.record{t}.cell(3) = struct('node',1,'record','v');
    neuron.pp{t}.IClamp = struct('node',1,'times',[-200 30,32.5],'amp', [hstep(t) hstep(t)+ostruct.cstep hstep(t)]); %n,del,dur,amp
    eucl{t} = eucl_tree(tree{t});
end
[out, ~] = t2n(neuron,tree,'-w-q-d');
tim = out.t;
tw = NaN(numel(tree),4,max(cellfun(@(x) numel(x{4}),CaNodes)));
maxcai = tw;


mCai = cell(numel(tree),numel(CaNodes{1}));
stdCai = mCai;
mV = mCai;
stdV = mCai;

for t = 1:numel(tree)
    plotcaivals{t} = NaN(numel(tree{t}.X),1);
    for x = 1:numel(nodes{t})
        plotcaivals{t}(nodes{t}(x)) = max(out.record{t}.cell.(ostruct.cai){nodes{t}(x)});
    end
    % as not all nodes were recorded, interpolate the value for the nodes
    % between the recorded ones (only affects plotting the tree, not the
    % data graphs)
    for x = 1:size(ipar{t},1)
        plotcaivals{t}(ipar{t}(x,:)) = interp1(plen{t}(intersect(nodes{t},ipar{t}(x,:))),plotcaivals{t}(intersect(nodes{t},ipar{t}(x,:))),plen{t}(ipar{t}(x,:)),'pchip');
    end
    for f =1:numel(CaNodes{t})
        for ff = 1:numel(CaNodes{t}{f})
            caivec = out.record{t}.cell.(ostruct.cai){CaNodes{t}{f}(ff)};
            if ~isempty(caivec)  % cai was existent at that node 
                opts = optimset('MaxFunEvals',50000,'MaxIter',10000);

                [maxcai(t,f,ff),ind(1)] = max(caivec);
                ind(2) = find(caivec(ind(1):end) <= caivec(end)+ 1e-7,1,'first')+ind(1)-1;
                if diff(ind) > 2
                    x = tim(ind(1):ind(2))-tim(ind(1));
                    yx = caivec(ind(1):ind(2))-caivec(ind(2));
                    y1 = @(b,x) b(1).*exp(-b(2).*x); % single exp function
                    y2 = @(b,x) b(1).*exp(-b(2).*x)+b(3).*exp(-b(4).*x);  % double exp function
                    B = fminsearch(@(b) sum((y2(b,x) - yx).^2),rand(4,1),opts); % do double exp fit
                    if B(1) < 0 || B(3) < 0  % if double exp fit fails (neg amplitude)
                        B = fminsearch(@(b) sum((y1(b,x) - yx).^2),rand(2,1),opts); % do single exp fit
                        tw(t,f,ff) = 1/B(2);
                    else
                        tw(t,f,ff) = (B(1)/B(2)+B(3)/B(4))/(B(1)+B(3)); % amplitude weighted tau as in stocca 2008
                    end
                else
                    tw(t,f,ff) = NaN;
                end
            else
                maxcai(t,f,ff) = NaN;
                tw(t,f,ff) = NaN;
            end
        end
        mCai{t,f} = mean(cat(2,out.record{t}.cell.(ostruct.cai){CaNodes{t}{f}})*1e6,2)';
        stdCai{t,f} = std (cat(2,out.record{t}.cell.(ostruct.cai){CaNodes{t}{f}})*1e6,[],2)';
        mV{t,f} = mean(cat(2,out.record{t}.cell.v{CaNodes{t}{f}}),2)';
        stdV{t,f} = std (cat(2,out.record{t}.cell.v{CaNodes{t}{f}}),[],2)';
    end
end

save(fullfile(targetfolder_data,sprintf('Exp_CaDyn_%s.mat',neuron.experiment)),'plotcaivals','nodes','neuron','tree','CaNodes','mCai','stdCai','mV','stdV','tim','tw','maxcai')
