session = {'20171214', '20171215', '20171218', '20171219'};
%1: 1st session with the cue at 0 or 200ms
%2: 1st session with the cue at 0 or 100ms
%3: 1st session with the cue at 0 or 50ms fp at 550ms
%4:  session with the cue at 0, 50 or 100ms and fp green at 600ms
Sessions = 3;
s=Sessions(1);
block =[1 2 1 1];
load(['D:\PostDoc\NeuronalData\Iconic\Monty\logs\Multicurve_Monty' session{s} '_B' num2str(block(s)) '.mat'])

load(['D:\PostDoc\NeuronalData\Iconic\Monty\extractdata\' session{s} '\EVT_Monty' session{s} '_Block-' num2str(block(s)) '.mat'])
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
CueT = LOG.NoCuedur(trials);
Poss_CueT = unique(CueT);
is1long = ones(numel(trials),1);%LOG.conn_codes(trials,1)-1;%

good_chan = [1:10 12:24 27:30 33:42];
avg_targ = zeros(numel(Poss_CueT),1373);
avg_nontarg = zeros(numel(Poss_CueT),1373);
avg_distr = zeros(6,1373);
for array = 1:2
    f1 = figure;
    f2 = figure;
    for n = (array-1)*24+1:array*24
        for s = Sessions
            load(['D:\PostDoc\NeuronalData\Iconic\Monty\extractdata\' session{s} '\Xtract_Monty' session{s} '_Block-' num2str(block(s)) '_' num2str(n) '.mat'])
            e = Env{1};
            e = e';
            figure(f1);subplot(4,6,n-(array-1)*24);hold on
            for ct = 1:numel(Poss_CueT)
                target = nanmean(e(RFonTargCurve & CueT==Poss_CueT(ct),:));
                nontarget = nanmean(e(RFonTargCurve'==0 & CueT'==Poss_CueT(ct) & is1long ,:));
                plot(tb,smooth(target,20)); plot(tb,smooth(nontarget,20),'r')
                xlim([-.2 1])
                if ~isempty(find(good_chan==n))
                    avg_targ(ct,:) = avg_targ(ct,:)+target;
                    avg_nontarg(ct,:) = avg_nontarg(ct,:)+nontarget;
                end
            end
            figure(f2);subplot(4,6,n-(array-1)*24);hold on
            for targn = 1:6
                data = nanmean(e(TargNum'==targn & is1long & CueT'==Poss_CueT(2),:));
                plot(tb, smooth(data,20)); xlim([-.2 1])
                if ~isempty(find(good_chan==n))
                    avg_distr(targn,:) = avg_distr(targn,:)+data;
                end
            end
            
            
            
        end
    end
end

figure;hold on;
for i = 1:numel(Poss_CueT)
    subplot(1, numel(Poss_CueT), i); hold on
    plot(tb,smooth(avg_targ(i,:),20), 'Color', [0 0 1]); plot(tb,smooth(avg_nontarg(i,:),20), 'Color', [1 0 0]); xlim([-.2 1])
    line([1e-3*Poss_CueT(i) 1e-3*Poss_CueT(i)], [-1 1], 'Color', 'k')
    line([.6 .6], [-1 1], 'Color', 'k')
    ylim([min(avg_targ(i,:)) max(avg_targ(i,:))])
end
% plot(tb,smooth(avg_targ(3,:),20), 'Color', [0 0 0.6]); plot(tb,smooth(avg_nontarg(3,:),20), 'Color', [0.6 0 0]); xlim([-.2 1])

% figure; plot(tb,smooth(avg_targ-avg_nontarg,20), 'g'); xlim([-.2 1])
figure; hold on
for targn = 1:6
    plot(tb, smooth(avg_distr(targn,:),20)); xlim([-.2 1])
end