function aGC_Ca_proportions_Elliot(neuron,tree)

elecnode = 1;

neuron.params.prerun = 300;
neuron.params.skiprun = 0;
neuron.params.tstop = 100;
neuron.params.dt= 1;  % 1 ms time step ist vollkommen ausreichend
neuron.params.cvode = 0;


for t = 1:numel(tree)
%     neuron.mech{t}.all.Cav12.kf = 1e9; %inactivate L-type inactivation
%     neuron.mech{t}.all.Cav13.kf = 1e9; %inactivate L-type inactivation
    fields = fieldnames(neuron.mech{t});
    for f = 1:numel(fields)
        if isfield(neuron.mech{t}.(fields{f}),'Cav12')
            neuron.mech{t}.(fields{f}).Cav12.kf = 1e9;%inactivate L-type inactivation
        end
        if isfield(neuron.mech{t}.(fields{f}),'Cav13')
            neuron.mech{t}.(fields{f}).Cav13.kf = 1e9;%inactivate L-type inactivation
        end
    end
    neuron.mech{t}.all.ca_ion.cao0 = 5; %extracellular barium concentration (calcium substitute)
    neuron.pp{t}.SEClamp = struct('node',elecnode,'times',[-200 20 20+60],'amp', [-100,0,-100],'rs',0.4);  % resistance was compensated by 80% CHECK IMPACT OF ELECTRODE RESISTANCE
    %         nneuron{s}.record{t}.cell = struct('record','v','node',elecnode);
    neuron.record{t}.SEClamp = struct('record','i','node',elecnode);
    OMLnode(t) = find(tree{t}.R==find(strcmp(tree{t}.rnames,'adendOML')),1);
    neuron.record{t}.cell = struct('node',[elecnode,OMLnode],'record',{'cai','v'});
end

% nneuron{1} = t2n_blockchannel(neuron,{'except','pas','Cav32','Ca2'},100);
% nneuron{2} = t2n_blockchannel(neuron,'Cav32',100);
% nneuron{3} = t2n_blockchannel(nneuron{1},'Ca2',100,[],{'gncabar'});
% nneuron{4} = t2n_blockchannel(nneuron{1},'Ca2',100,[],{'glcabar'});

nneuron{1} = t2n_blockchannel(neuron,{'except','pas','Cav22','Cav32','Cav12','Cav13'},100);
nneuron{2} = t2n_blockchannel(nneuron{1},{'Cav32'},100);
nneuron{3} = t2n_blockchannel(nneuron{1},'Cav22',100);
nneuron{4} = t2n_blockchannel(nneuron{1},{'Cav12','Cav13'},100);
nneuron{5} = t2n_blockchannel(nneuron{1},{'Cav32','Cav12','Cav13','Cav22'},100);

[out, ~] = t2n(nneuron,tree,'-q-d-w');
if any(cellfun(@(x) x.error,out))
    return
end

tvec = out{1}.t;

cai = zeros(numel(tvec),numel(tree)); curr_Nblock = cai; curr_Lblock = cai; gak = cai; gabk = cai;gakbar = cai; gabkbar = cai;v = cai;

for t = 1:numel(tree)
%         figure; hold all
%         plot(out{1}.t,out{1}.record{t}.SEClamp.i{elecnode})
%         plot(out{2}.t,out{2}.record{t}.SEClamp.i{elecnode})
%         plot(out{3}.t,out{3}.record{t}.SEClamp.i{elecnode})
%         plot(out{4}.t,out{4}.record{t}.SEClamp.i{elecnode})
%         plot(out{5}.t,out{5}.record{t}.SEClamp.i{elecnode})
%     figure, hold all
%     plot(out{1}.t,out{1}.record{t}.cell.g_Cav22{elecnode})
%     plot(out{1}.t,out{1}.record{t}.cell.glca_Ca2{elecnode})
%     plot(out{1}.t,out{1}.record{t}.cell.g_Cav32{elecnode})
    curr_all(:,t) = out{1}.record{t}.SEClamp.i{elecnode} - out{5}.record{t}.SEClamp.i{elecnode};
    curr_L(:,t) = out{1}.record{t}.SEClamp.i{elecnode} - out{4}.record{t}.SEClamp.i{elecnode};
    curr_N(:,t) = out{1}.record{t}.SEClamp.i{elecnode} - out{3}.record{t}.SEClamp.i{elecnode};
    curr_T(:,t) = out{1}.record{t}.SEClamp.i{elecnode} - out{2}.record{t}.SEClamp.i{elecnode};
    
    %     cai(:,t) = out{1}.record{t}.cell.cai{elecnode};
    %     gabk(:,t) = out{1}.record{t}.cell.gabk_BK{elecnode};
    %     gak(:,t) = out{1}.record{t}.cell.gak_BK{elecnode};
    %     gabkbar(:,t) = out{1}.record{t}.cell.gabkbar_BK{elecnode};
    %     gakbar(:,t) = out{1}.record{t}.cell.gakbar_BK{elecnode};
    %     v(:,t) = out{1}.record{t}.cell.v{elecnode};
%     gnca(:,t) = out{1}.record{t}.cell.g_Cav22{elecnode};
%     glca(:,t) = out{1}.record{t}.cell.glca_Ca2{elecnode};
%     gtca(:,t) = out{1}.record{t}.cell.g_Cav32{elecnode};
    
%     gnca_OML(:,t) = out{1}.record{t}.cell.g_Cav22{OMLnode(t)};
%     glca_OML(:,t) = out{1}.record{t}.cell.glca_Ca2{OMLnode(t)};
%     gtca_OML(:,t) = out{1}.record{t}.cell.g_Cav32{OMLnode(t)};
end

ind = out{1}.t > 30 & out{1}.t < 80 ;
% max(curr_L(ind,:)./curr_all(ind,:),[],1)*100
% max(curr_N(ind,:)./curr_all(ind,:),[],1)*100
% max(curr_T(ind,:)./curr_all(ind,:),[],1)*100

% figure;hold all,
% p1= plot(tvec,gnca./(gnca+gtca+glca),'g');
% ylabel('X-type rel. conductance contribution to Ca')
% xlabel('time [ms]')
% p2= plot(tvec,glca./(gnca+gtca+glca),'b');
% p3 = plot(tvec,gtca./(gnca+gtca+glca),'k');
% legend([p1(1),p2(1),p3(1)],'N-type','L-type','T-type')
% 
% figure;hold all,
% p1= plot(tvec,gnca,'g');
% ylabel('X-type conductance')
% xlabel('time [ms]')
% p2= plot(tvec,glca,'b');
% p3 = plot(tvec,gtca,'k');
% p4 = plot(tvec,gnca+gtca+glca,'c--');
% legend([p1(1),p2(1),p3(1),p4(1)],'N-type','L-type','T-type','Total')

figure;hold all,
p1= plot(tvec,curr_N./curr_all,'g');
ylabel('X-type rel. current contribution to Ca')
xlabel('time [ms]')
p2= plot(tvec,curr_L./curr_all,'b');
p3 = plot(tvec,curr_T./curr_all,'k');
legend([p1(1),p2(1),p3(1)],'N-type','L-type','T-type')

figure;hold all,
p1= plot(tvec,curr_N,'g');
ylabel('X-type current [nA]')
xlabel('time [ms]')
p2= plot(tvec,curr_L,'b');
p3 = plot(tvec,curr_T,'k');
p4 = plot(tvec,curr_all,'c--');
legend([p1(1),p2(1),p3(1),p4(1)],'N-type','L-type','T-type','Total')

% figure;hold all,
% p1= plot(tvec,gnca_OML,'g');
% ylabel('X-type conductance in OML')
% xlabel('time [ms]')
% p2= plot(tvec,glca_OML,'b');
% p3 = plot(tvec,gtca_OML,'k');
% p4 = plot(tvec,gnca_OML+gtca_OML+glca_OML,'c--');
% legend([p1(1),p2(1),p3(1),p4(1)],'N-type','L-type','T-type','Total')

% 
% display(sprintf('Iberitoxin-sensitive current: %g +- %g nA (s.e.m.)',mean(abscurr_ibersens),std(abscurr_ibersens)/sqrt(numel(tree))));
% display(sprintf('Paxilline-sensitive current: %g +- %g nA (s.e.m.)',mean(abscurr_paxisens),std(abscurr_paxisens)/sqrt(numel(tree))));
% display(sprintf('Part of typeI BK in total BK current: %g +- %g %% (s.e.m.)',mean(ratio_type1),std(ratio_type1)/sqrt(numel(tree))));
