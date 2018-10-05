session = {'20171214', '20171215', '20171218', '20171219', '20180109', '20180110'};
%1:  session  with the cue at 0 or 200ms
%2:  session  with the cue at 0 or 100ms
%3:  session  with the cue at 0 or 50ms and           fp green at 550ms
%4:  session  with the cue at 0, 50 or 100ms and      fp green at 600ms
%5-6:sessions with the cue at 0, 100 or 200ms and fp green at 900ms
ana_sessions = 5:6;
block ={1, 2, 1, 1, 1, 1:3};
load(['D:\PostDoc\NeuronalData\Iconic\Monty\logs\Multicurve_Monty' session{ana_sessions(1)} '_B' num2str(block{ana_sessions(1)}(1)) '.mat'])

load(['D:\PostDoc\NeuronalData\Iconic\Monty\extractdata\' session{ana_sessions(1)} '\EVT_Monty' session{ana_sessions(1)} '_Block-' num2str(block{ana_sessions(1)}(1)) '.mat'])
SF = EVENT.strms(1).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);


close all

trials = find(LOG.target_presented);
CueT = LOG.NoCuedur(trials);
Poss_CueT = unique(CueT);
% is1long = ones(numel(trials),1);%LOG.conn_codes(trials,1)-1;%

% good_chan = [1:8 10:24 27:30 33:42];
avg_targ = zeros(numel(Poss_CueT),1373);
avg_nontarg = zeros(numel(Poss_CueT),1373);
avg_distr = zeros(6,numel(Poss_CueT),1373);
for array = 1:2
    f1 = figure;
    f2 = figure;
    for n = (array-1)*24+1:array*24
        CondTarget = cell(numel(Poss_CueT),6);
        for s = ana_sessions
            for b = block{s}
                load(['D:\PostDoc\NeuronalData\Iconic\Monty\extractdata\' session{s} '\Xtract_Monty' session{s} '_Block-' num2str(block{s}(b)) '_' num2str(n) '.mat'])
                e = Env{1};
                e = e';
                load(['D:\PostDoc\NeuronalData\Iconic\Monty\logs\Multicurve_Monty' session{s} '_B' num2str(block{s}(b)) '.mat'])
                trials = find(LOG.target_presented);
                if numel(trials)>size(e,1)
                    trials = trials(1:size(e,1));
                    disp('Warning: More trials in the Log file than in the data')
                end
                TargNum = LOG.targ_num(trials);
                RFonTargCurve = TargNum==1;
                CueT = LOG.NoCuedur(trials);
                Poss_CueT = unique(CueT);
                is1long = ones(numel(trials),1);
                
                
                for ct = 1:numel(Poss_CueT)
                    for targn = 1:6
                        TargTrials = e(TargNum'==targn & is1long & CueT'==Poss_CueT(ct),:);
                        if isempty(CondTarget{ct,targn})
                            CondTarget{ct,targn} = TargTrials;
                        else
                            CondTarget{ct,targn} = [CondTarget{ct,targn};TargTrials];
                        end
                    end
                end
            end
        end
        target = nanmean(CondTarget{1,1});
        SNR = (max(target(tb>0.05 & tb<0.1))-mean(target(tb>-.2 & tb<0)))./mean(std(e(RFonTargCurve & CueT==Poss_CueT(ct),tb>-.2 & tb<0),1,2));
        isgood = SNR>.5;
        if isgood
            figure('Position',[100,100,1400,700]);
            for ct = 1:numel(Poss_CueT)
                target = nanmean(CondTarget{ct,1});
                nontarget = .2*(nanmean(CondTarget{ct,2})+nanmean(CondTarget{ct,3})+nanmean(CondTarget{ct,4})+nanmean(CondTarget{ct,5})+nanmean(CondTarget{ct,6}));
                avg_targ(ct,:) = avg_targ(ct,:)+target;
                avg_nontarg(ct,:) = avg_nontarg(ct,:)+nontarget;
                for targn = 1:6
                    avg_distr(targn,ct,:) = avg_distr(targn,ct,:)+reshape(nanmean(CondTarget{ct,targn}),1, 1, numel(nanmean(CondTarget{ct,targn})));
                end
                subplot(2,numel(Poss_CueT),ct);hold on
                plot(tb,smooth(target,20), 'Color', 'b'); plot(tb,smooth(nontarget,20),'Color', 'r')
                xlim([-.2 1]) 
                if ct==1
                    legend('target', 'distractor')
                    title(['Chan' num2str(n)])
                end
                subplot(2,numel(Poss_CueT),ct+numel(Poss_CueT));hold on
                for targn = 1:6
                    plot(tb,smooth(nanmean(CondTarget{ct,targn}),20))
                end
                xlim([-.2 1]) 
                title(['Cue at ' num2str(Poss_CueT(ct)) 'ms'])
                if ct==1
                    legend('target', 'same pair distractor')
                end
            end
            pause(0.05)
        end
        
    end
    
end
%%             
figure('Position',[100,100,1400,400]);hold on;
for i = 1:numel(Poss_CueT)
    subplot(1, numel(Poss_CueT), i); hold on
    plot(tb,smooth(avg_targ(i,:),20), 'Color', [0 0 1]); plot(tb,smooth(avg_nontarg(i,:),20), 'Color', [1 0 0]); xlim([-.2 1])
    line([1e-3*Poss_CueT(i) 1e-3*Poss_CueT(i)], [-1 1], 'Color', 'k')
    line([1e-3*Poss_CueT(i)+.15 1e-3*Poss_CueT(i)+.15], [-1 1], 'Color', 'g')
    line([1e-3*Poss_CueT(i)+.35 1e-3*Poss_CueT(i)+.35], [-1 1], 'Color', 'y')
    ylim([min(avg_targ(1,:)) max(avg_targ(1,:))+1e-5])
    title(['Cue at ' num2str(Poss_CueT(i)) 'ms'])
    if i==1
        legend('target', 'distractor')
    end
end
% plot(tb,smooth(avg_targ(3,:),20), 'Color', [0 0 0.6]); plot(tb,smooth(avg_nontarg(3,:),20), 'Color', [0.6 0 0]); xlim([-.2 1])

% figure; plot(tb,smooth(avg_targ-avg_nontarg,20), 'g'); xlim([-.2 1])
figure('Position',[100,100,1400,400]); hold on
for ct = 1:numel(Poss_CueT)
    subplot(1, numel(Poss_CueT), ct); hold on
    for targn = 1:6
        plot(tb, smooth(avg_distr(targn,ct,:),20)); xlim([-.2 1])
    end
    line([1e-3*Poss_CueT(ct) 1e-3*Poss_CueT(ct)], [-1 1], 'Color', 'k')
    line([1e-3*Poss_CueT(ct)+.1  1e-3*Poss_CueT(ct)+.1 ], [-1 1], 'Color', 'g')
    line([1e-3*Poss_CueT(ct)+.35 1e-3*Poss_CueT(ct)+.35], [-1 1], 'Color', 'y')
    ylim([min(avg_targ(1,:)) max(avg_targ(1,:))])
    title(['Cue at ' num2str(Poss_CueT(ct)) 'ms'])
    if ct==1
        legend('target', 'same pair distractor')
    end
end