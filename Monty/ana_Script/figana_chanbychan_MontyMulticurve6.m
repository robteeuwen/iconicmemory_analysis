session = {'20171204', '20171205', '20171206', '20171207', '20171208', '20171211', '20171212'};
%3: session with 6 pairs
%4:5: session with targets closer together 
%6: 1st session with 1short and 1 long segment per pair
%7: 1st session with initial segment
s = 7;
block =2;
load(['D:\PostDoc\NeuronalData\Iconic\Monty\logs\Multicurve_Monty' session{s} '_B' num2str(block) '.mat'])

load(['D:\PostDoc\NeuronalData\Iconic\Monty\extractdata\' session{s} '\EVT_Monty' session{s} '_Block-' num2str(block) '.mat'])
SF = EVENT.strms(1).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);




trials = find(LOG.target_presented);
TargNum = LOG.targ_num(trials);
% sideX = LOG.RFx(trials,:)>0;
% sideY = LOG.RFy(trials,:)<0;
% isinRF = (sideX+sideY)==2;
% [targinRF,b] = find(isinRF');
RFonTargCurve = TargNum==1;
is1long = LOG.conn_codes(trials,1)-1;

good_chan = [1:10 12:24 27:30 33:42];
avg_targ = zeros(1,1373);
avg_nontarg = zeros(1,1373);
avg_distr = zeros(6,1373);
for array = 1:2
    f1 = figure; 
    f2 = figure;
    for n = (array-1)*24+1:array*24
        load(['D:\PostDoc\NeuronalData\Iconic\Monty\extractdata\' session{s} '\Xtract_Monty' session{s} '_Block-' num2str(block) '_' num2str(n) '.mat'])
        e = Env{1};
        e = e';
        target = nanmean(e(RFonTargCurve,:));
        nontarget = nanmean(e(RFonTargCurve'==0 & is1long,:));
        figure(f1);subplot(4,6,n-(array-1)*24);plot(tb,smooth(connected,20));hold on; plot(tb,smooth(unconnected,20),'r')
        xlim([-.2 1])
        figure(f2);subplot(4,6,n-(array-1)*24);hold on
        for targn = 1:6
            data = nanmean(e(TargNum'==targn & is1long ,:));
            plot(tb, smooth(data,20)); xlim([-.2 1])
        end
        if ~isempty(find(good_chan==n))
            avg_targ = avg_targ+target;
            avg_nontarg = avg_nontarg+nontarget;
            for d = 1:6
                avg_distr(d,:) = avg_distr(d,:)+nanmean(e(TargNum'==d  & is1long,:));
            end
        end
    end
    
end

figure; plot(tb,smooth(avg_targ,20));hold on; plot(tb,smooth(avg_nontarg,20),'r'); xlim([-.2 1])
figure; plot(tb,smooth(avg_targ-avg_nontarg,20), 'g'); xlim([-.2 1])
figure; hold on 
for targn = 1:6
    plot(tb, smooth(avg_distr(targn,:),20)); xlim([-.2 1])
end