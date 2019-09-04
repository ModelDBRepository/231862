function aGC_Shruti_current(neuron,tree,ostruct)
% if nargin < 5 || ~isfield(ostruct,'holding_voltage')
ostruct.holding_voltage = -80;
% end
elecnode = 1;

neuron.params.prerun = 300;
neuron.params.skiprun = 0;
neuron.params.tstop = 920;
neuron.params.dt= 1;  % 1 ms time step ist vollkommen ausreichend
neuron.params.cvode = 0;


for t = 1:numel(tree)
    neuron.pp{t}.SEClamp = struct('node',elecnode,'times',[-500 50 50+40 50+40+800],'amp', [-80,-40,+100,-80],'rs',6);  % CHECK IMPACT OF ELECTRODE RESISTANCE
    %         nneuron{s}.record{t}.cell = struct('record','v','node',elecnode);
    neuron.record{t}.SEClamp = struct('record','i','node',elecnode);
    neuron.record{t}.cell = struct('node',1,'record',{'gak_BK','gabk_BK','gakbar_BK','gabkbar_BK','cai','v'});
end

nneuron{1} = t2n_blockchannel(neuron,{'na8st'},100);
nneuron{2} = t2n_blockchannel(neuron,{'na8st','BK'},100,[],{'gbar','gakbar'});   % block with iberitoxin
nneuron{3} = t2n_blockchannel(neuron,{'na8st','BK'},100);         % block with paxilline

[out, ~] = t2n(nneuron,tree,'-q-d-w');
if any(cellfun(@(x) x.error,out))
    return
end

tvec = out{1}.t;

cai = zeros(numel(tvec),numel(tree)); curr_paxisens = cai; curr_ibersens = cai; gak = cai; gabk = cai;gakbar = cai; gabkbar = cai;v = cai;

for t = 1:numel(tree)
    curr_ibersens(:,t) = out{1}.record{t}.SEClamp.i{elecnode} - out{2}.record{t}.SEClamp.i{elecnode};
    curr_paxisens(:,t) = out{1}.record{t}.SEClamp.i{elecnode} - out{3}.record{t}.SEClamp.i{elecnode};
    cai(:,t) = out{1}.record{t}.cell.cai{elecnode};
    gabk(:,t) = out{1}.record{t}.cell.gabk_BK{elecnode};
    gak(:,t) = out{1}.record{t}.cell.gak_BK{elecnode};
    gabkbar(:,t) = out{1}.record{t}.cell.gabkbar_BK{elecnode};
    gakbar(:,t) = out{1}.record{t}.cell.gakbar_BK{elecnode};
    v(:,t) = out{1}.record{t}.cell.v{elecnode};
end

abscurr_ibersens = max(curr_ibersens(tvec>590&tvec<890,:),[],1);
abscurr_paxisens = max(curr_paxisens(tvec>590&tvec<890,:),[],1);
ratio_type1 = abscurr_ibersens ./ abscurr_paxisens * 100;

figure;hold on
% plot(tvec,out{1}.record{t}.SEClamp.i{elecnode},'b')
plot(tvec,curr_ibersens,'k')
plot(tvec,curr_paxisens,'Color',[0.5 0.5 0.5])
xlabel('Time [ms]')
ylabel('Subtracted current [nA]')
legend('Iberitoxin-sensitive','Paxilline-sensitive')
xlim([80 920])

figure
plot(tvec,cai)
xlabel('Time [ms]')
ylabel('Calcium concentration [mM?]')
xlim([80 920])

figure
plot(tvec,v)
xlabel('Time [ms]')
ylabel('Membrane potential [mV]')
xlim([80 920])

figure;hold on
p1=plot(tvec,gak,'k');
p2=plot(tvec,gabk,'Color',[0.5 0.5 0.5]);
plot(tvec,gakbar,'r')
plot(tvec,gabkbar,'Color',[1 0.2 0.2])
xlabel('Time [ms]')
ylabel('Conductance [S/cm²]')
legend([p1(1) p2(1)],'Type I BK','Type II BK')
xlim([80 920])

fprintf('Iberitoxin-sensitive current: %g +- %g nA (s.e.m.)\n',mean(abscurr_ibersens),std(abscurr_ibersens)/sqrt(numel(tree)));
fprintf('Paxilline-sensitive current: %g +- %g nA (s.e.m.)\n',mean(abscurr_paxisens),std(abscurr_paxisens)/sqrt(numel(tree)));
fprintf('Part of typeI BK in total BK current: %g +- %g %% (s.e.m.)\n',mean(ratio_type1),std(ratio_type1)/sqrt(numel(tree)));