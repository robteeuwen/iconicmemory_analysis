% analysis function for the reboot of the cueing part of the iconic memory
% project
% sessions 94 - 122
% written by Catherine Wacongne (catherine.waco@gmail.com)


close all
clearvars
dbstop if error
plot_stim_dur =1;
plotask = 1;
plot_chann= 1;
plot_corr = 0;

integ_dur = 0.30;
smoothness = 30;
baseline = [-.5 -.3];
x_limits = [-.5 1.2];

[extractdir, sessions, rawdir, SaveDir] = infoDir_MontyIcoMemo;
info = Log_MontyIcoMemo;
load([extractdir, sessions{end} '\EVT_', info(numel(sessions)).Tankname,'_block-',num2str(info(numel(sessions)).goodblocs(1))])

SF = EVENT.strms(1).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);
StimCondStr = {'Short', 'Long'};
Cond = cell(96,2,2,2);
ci = zeros(96,2);


%94-97 reboot of the cueing paradigm with the constant 3 curves
%98-105  adds delays 50 and 100ms after the end of stim
%106-114 adds delays 160ms after the end of stim, longer trials + longer
%115- 122 disapearing segm in RF
Sessions = 106:114;%
SNR = zeros(1,48);

for array = 2
    ana_chan = (array-1)*24+1:array*24;
    if array==2
        ana_chan = 26:45;
    end
    for stimcond=1
        aregood = zeros(max(Sessions),48);
        for n = ana_chan %48
            
            % Extract the trials from all blocs and conditions for each channel
            
            AllTrials = [];
            AllCond = [];
            days =[];
            for s = Sessions%60;%28:29%10:23%10:21;%10:12%6:numel(sessions)
                
                for bloc = 1:numel(info(s).goodblocs)
                    %Read in the extracted data form the mapped channel
                    load([extractdir sessions{s},'\Xtract_',info(s).Tankname,'_Block-',num2str(info(s).goodblocs(bloc)),'_',num2str(n)]);
                    load([extractdir, sessions{s} '\EVT_', info(s).Tankname,'_block-',num2str(info(s).goodblocs(bloc))])

                    e = Env{1};
                    load([rawdir,'ICOmemo_' info(s).Tankname '_B',num2str(info(s).goodblocs(bloc))])
                    trials =find(LOG.target_presented>0);
                    if numel(trials)>size(e,2)
%                         keyboard
                        trials = trials(1:size(e,2));
                        
                    end
                    days = [days, s*ones(1,size(e,2))];
                    if isempty(AllTrials)
                        AllTrials=e;
                        %                         connected = LOG.conn_codes(trials,1);
                        targ = LOG.targ_num(trials);
                        correct = LOG.correct(trials)+1;
                        cue_time = LOG.CueDelay(trials);
                        stim_dur = LOG.Stim_dur(trials); % Stimulus duration
                        islong = LOG.conn_codes(trials,1)'-1;
                        angle = LOG.angle(trials);
                    else
                        AllTrials = [AllTrials e];
                        %                         connected = [connected LOG.conn_codes(trials,1)];
                        targ = [targ LOG.targ_num(trials)];
                        correct = [correct LOG.correct(trials)+1];
                        cue_time = [cue_time LOG.CueDelay(trials)];
                        stim_dur = [stim_dur LOG.Stim_dur(trials)]; % Stimulus duration
                        islong = [islong LOG.conn_codes(trials,1)'-1];
                        angle = [angle LOG.angle(trials)];
                    end
                    
                end
            end
            condcheck = 0*targ;
            Poss_dur = unique(stim_dur);
            fix_stim_dur = Poss_dur(1)*1.e-3;
            Poss_Cue_Time = unique(cue_time);
            Poss_Angle = unique(angle);      
            bf = find(tb<baseline(2));
            if s<115
                pf = find(tb > 0.05-LOG.Precue_dur(1)/1e3 & tb < 0.1-LOG.Precue_dur(1)/1e3);
                pf2 = find(tb > 0.05-LOG.Precue_dur(1)/1e3 & tb < 0.09-LOG.Precue_dur(1)/1e3);
            else
                pf = find(tb > 0.05 & tb < 0.1);
                pf2 = find(tb > 0.05 & tb < 0.09);
            end
                
            hitNorm =(targ==1 & correct==2);
            %      AllTrials = NormalizeMUA(AllTrials, n, days,hitNorm ,bf, pf);
            [AllTrials2,eraseday, isgood] = NormalizeActiv(AllTrials, n, days,hitNorm ,bf, pf);
            %     SNR(n) = nanmean(reshape(AllTrials2(pf2,hitNorm==1),1,[]))/nanmean(nanstd(AllTrials2(bf,:)));
            aregood(:,n) = isgood;
            SNR(n) = nanmax(nanmean(AllTrials2(pf2,hitNorm==1),2))/nanmean(nanstd(AllTrials2(bf,:),0,2));
            %     if mean(eraseday)>.5
            %         SNR(n) = 0;
            %     end
            %     figure; subplot(1,2,1);imagesc(AllTrials);subplot(1,2,2); imagesc(AllTrials2);pause(0.1)
            AllTrials = AllTrials2;
            %      Fix_sd = 60;
            ind_a = cell(1,numel(Poss_Angle));
            for h = 1:2
                for ct = 1:numel(Poss_Cue_Time)
                    
                    % average per condition: target tree - targ
                    for ang = 1:numel(Poss_Angle)
                        if plot_corr
                            ind_a{ang} = find(targ==1  & correct==h & cue_time==Poss_Cue_Time(ct) & angle==Poss_Angle(ang));%
                        else
                            ind_a{ang} = find(targ==1  & cue_time==Poss_Cue_Time(ct) & angle==Poss_Angle(ang));
                        end
                    end
                    mintrialspercond = min(cellfun(@numel,ind_a));
                    ind = [];
                    for ang = 1:numel(Poss_Angle)
                        rd = randperm(numel(ind_a{ang}));
                        ind = [ind ind_a{ang}(rd(1:mintrialspercond))];
                    end
                    ttree_t{n,h,ct} = nanmean(AllTrials(:,ind),2);
                    %                     if ~isempty(find(condcheck(ind))); keyboard;end
                    %                     condcheck(ind) = 2;
                    
                    
                    % targ tree - distr
                    for ang = 1:numel(Poss_Angle)
                        if plot_corr
                            ind_a{ang} = find(targ==2 & correct==h  & cue_time==Poss_Cue_Time(ct)& angle==Poss_Angle(ang));
                        else
                            ind_a{ang} = find(targ==2  & cue_time==Poss_Cue_Time(ct)& angle==Poss_Angle(ang));
                        end
                    end
                    mintrialspercond = min(cellfun(@numel,ind_a));
                    ind = [];
                    for ang = 1:numel(Poss_Angle)
                        rd = randperm(numel(ind_a{ang}));
                        ind = [ind ind_a{ang}(rd(1:mintrialspercond))];
                    end
                    ttree_d{n,h,ct} = nanmean(AllTrials(:,ind),2);
                    %                     if ~isempty(find(condcheck(ind))); keyboard;end
                    %                     condcheck(ind) = 3;
                    
                    
                    
                    % distr tree - connected
                    for ang = 1:numel(Poss_Angle)
                        if plot_corr
                            ind_a{ang} = find(islong==1 & targ~=1  & correct==h & cue_time==Poss_Cue_Time(ct)& angle==Poss_Angle(ang));
                        else
                            ind_a{ang} = find(islong==1  & targ~=1 & cue_time==Poss_Cue_Time(ct)& angle==Poss_Angle(ang));
                        end
                    end
                    mintrialspercond = min(cellfun(@numel,ind_a));
                    ind = [];
                    for ang = 1:numel(Poss_Angle)
                        rd = randperm(numel(ind_a{ang}));
                        ind = [ind ind_a{ang}(rd(1:mintrialspercond))];
                    end
                    dtree_t{n,h,ct} = nanmean(AllTrials(:,ind),2);
                    %                     if ~isempty(find(condcheck(ind))); keyboard;end
                    %                     condcheck(ind) = 4;
                    
                    % distr tree - unconnected
                    for ang = 1:numel(Poss_Angle)
                        if plot_corr
                            ind_a{ang} = find(islong==0 & targ~=2 & correct==h & cue_time==Poss_Cue_Time(ct)& angle==Poss_Angle(ang));
                        else
                            ind_a{ang} = find(islong==0 & targ~=2 & cue_time==Poss_Cue_Time(ct)& angle==Poss_Angle(ang));
                        end
                    end
                    mintrialspercond = min(cellfun(@numel,ind_a));
                    ind = [];
                    for ang = 1:numel(Poss_Angle)
                        rd = randperm(numel(ind_a{ang}));
                        ind = [ind ind_a{ang}(rd(1:mintrialspercond))];
                    end
                    dtree_d{n,h,ct} = nanmean(AllTrials(:,ind),2);
                    %                     if ~isempty(find(condcheck(ind))); keyboard;end
                    %                     condcheck(ind) = 5;
                    
                    ind = find( correct==h  & cue_time==Poss_Cue_Time(ct));
                    perf(h, ct) = numel(ind);
                    
                end
                %                 [hyp, ismodul(n),ci(n,:),stats] = ttest2(nanmean(AllTrials(pf2,ind1)),nanmean(AllTrials(pf2,ind2))) ;
                
                
            end
            
        end
        
        %% plot the conditions for each channel
        
        close all
        if plot_chann
            %     f1 =figure('Position', [100 50 1400 900]);
            f2 =figure('Position', [100 50 1400 900]);
            %n,h,sd,m
            
            AUC = zeros(min(numel(Poss_Cue_Time),6),48);
            for n = ana_chan
                
                figure(f2);clf
                for ct = 1:min(numel(Poss_Cue_Time),6)
                    
                    figure(f2)
                    subplot(2,4,ct);hold on; title([num2str(Poss_Cue_Time(ct)), ' ms, chan ' num2str(n)])
                    
                    
                    plot(tb, smooth(ttree_t{n,2,ct}, smoothness), 'b', 'LineWidth',1); try  ylim([min(smooth(ttree_t{n,2,ct},smoothness))-(max(smooth(ttree_t{n,2,ct},smoothness))-min(smooth(ttree_t{n,2,ct},smoothness))),  max(smooth(ttree_t{n,2,ct},smoothness))+.2*(max(smooth(ttree_t{n,2,ct},smoothness))-min(smooth(ttree_t{n,2,ct},smoothness)))]);end
                    plot(tb, smooth(ttree_d{n,2,ct}, smoothness), 'r', 'LineWidth',1);
                    plot(tb, smooth(dtree_d{n,2,ct}, smoothness), 'm', 'LineWidth',1);
                    [hAx,hLine1,hLine2] = plotyy(tb, smooth(dtree_t{n,2,ct}, smoothness), tb, [smooth(ttree_t{n,2,ct}-ttree_d{n,2,ct}, smoothness)';smooth(dtree_t{n,2,ct}-dtree_d{n,2,ct}, smoothness)']);
                    diffcu =smooth(ttree_t{n,2,ct}-dtree_t{n,2,ct},smoothness);
                    hAx(2).XLim = x_limits;
                    hAx(1).XLim = x_limits;
                    hLine1.LineWidth = 1;
                    hLine1.Color = [0 1 1];
                    try hAx(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
                    hLine2(1).Color=[0 1 0];
                    hLine2(2).Color=[1 .5 0];
                    hLine2(1).LineWidth = 1;
                    hLine2(2).LineWidth = 1;
                    
                    
                    axes(hAx(2));line(hAx(2).XLim,[0 0]);
                    if Poss_Cue_Time(ct)<0
                        p=patch([0.1 0.1+integ_dur 0.1+integ_dur 0.1],[-10 -10 20 20],'r');
                        t_relev = find(tb>(0.1) & tb<(0.1+integ_dur));
                    else
                        p=patch([fix_stim_dur+Poss_Cue_Time(ct)*1e-3, fix_stim_dur+Poss_Cue_Time(ct)*1e-3+integ_dur, fix_stim_dur+Poss_Cue_Time(ct)*1e-3+integ_dur, fix_stim_dur+Poss_Cue_Time(ct)*1e-3],[-10 -10 20 20],'r');
                        t_relev = find(tb>(fix_stim_dur+Poss_Cue_Time(ct)*1e-3) & tb<(fix_stim_dur+Poss_Cue_Time(ct)*1e-3 +integ_dur));
                    end
                    set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
                    diff =abs(smooth(ttree_t{n,2,ct}, smoothness)-smooth(ttree_d{n,2,ct}, smoothness)-smooth(dtree_t{n,2,ct}, smoothness)+smooth(dtree_d{n,2,ct}, smoothness));
                    diff = diff(t_relev);%diff(diff<0)=0;
                    %                     i = find(diff<0, 1);diff(i:end)=0;
                    AUC(ct,n) = sum(diff);
                    if Poss_Cue_Time(ct)>=0
                        line([fix_stim_dur+1e-3*Poss_Cue_Time(ct) fix_stim_dur+1e-3*Poss_Cue_Time(ct)], [-10 20])
                    end
                end
                
                subplot(2,4,8);hold on; plot(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC(:,n), 'r')
                xlabel('cue delay')
                pause(0.2)
                set(gcf,'PaperPositionMode','auto')
                print(gcf,'-dpng',[SaveDir 'ConnUnconn_Chan' num2str(n) '_s' num2str(days(1)) '_' num2str(days(end)) '.png'])
                
            end
        end
        %% plot performances
        good_chann = find(SNR>.5);
        
        figure;
        fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
        title('Behavioral performance on the iconic memory task');
        xlabel('Stimulus-cue interval (ms)')
        
        hold on;  plot(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)),  perf(2,1:min(numel(Poss_Cue_Time),6))./sum(perf(:,1:min(numel(Poss_Cue_Time),6))),'r')
        
        filename = ['performances_' StimCondStr{stimcond} 'Stim_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        %% plot the average over an array for all the conditions
        f2 = figure('Position', [100 50 1400 900]);f3 = figure('Position', [100 50 1400 900]);
        for ct = 1:min(numel(Poss_Cue_Time),6)
            % compute the average
            for h = 1:2
                ttree_t_avg{h,ct} = 0*ttree_t{ana_chan(1),1};n_ttree_t_avg(h,ct)=0;
                ttree_d_avg{h,ct} = 0*ttree_t{ana_chan(1),1};n_ttree_d_avg(h,ct)=0;
                dtree_t_avg{h,ct} = 0*ttree_t{ana_chan(1),1};n_dtree_t_avg(h,ct)=0;
                dtree_d_avg{h,ct} = 0*ttree_t{ana_chan(1),1};n_dtree_d_avg(h,ct)=0;
            end
            for n = good_chann
                for h = 1:2
                    try  ttree_t_avg{h,ct} =  nansum([ttree_t_avg{h,ct} ,ttree_t{n,h,ct}],2);end;if ~isnan(sum(ttree_t{n,h,ct}));n_ttree_t_avg(h,ct) = n_ttree_t_avg(h,ct)+1;end
                    try ttree_d_avg{h,ct} =  nansum([ttree_d_avg{h,ct} ,ttree_d{n,h,ct} ],2);end;if ~isnan(sum(ttree_d{n,h,ct}));n_ttree_d_avg(h,ct) = n_ttree_d_avg(h,ct)+1;end
                    
                    try dtree_t_avg{h,ct} =  nansum([dtree_t_avg{h,ct},dtree_t{n,h,ct}],2);end;if ~isnan(sum(dtree_t{n,h,ct}));n_dtree_t_avg(h,ct) = n_dtree_t_avg(h,ct)+1;end
                    try dtree_d_avg{h,ct} =  nansum([dtree_d_avg{h,ct},dtree_d{n,h,ct}],2);end;if ~isnan(sum(dtree_d{n,h,ct}));n_dtree_d_avg(h,ct) = n_dtree_d_avg(h,ct)+1;end
                end
                
            end
            for h = 1:2
                ttree_t_avg{h,ct} = ttree_t_avg{h,ct}/n_ttree_t_avg(h,ct);
                ttree_d_avg{h,ct} = ttree_d_avg{h,ct}/n_ttree_d_avg(h,ct);
                dtree_t_avg{h,ct} = dtree_t_avg{h,ct}/n_dtree_t_avg(h,ct);
                dtree_d_avg{h,ct} = dtree_d_avg{h,ct}/n_dtree_d_avg(h,ct);
            end
            
            % plot
            figure(f2);subplot(2,3,ct);hold on; title(['cued at ' num2str(Poss_Cue_Time(ct)), ' ms'])%
            plot(tb,smooth(ttree_t_avg{1,ct},smoothness),'b');
            plot(tb,smooth(ttree_d_avg{1,ct},smoothness),'r');
            plot(tb,smooth(dtree_t_avg{1,ct},smoothness),'c');
            plot(tb,smooth(dtree_d_avg{1,ct},smoothness),'m');
            if ct == 1; legend({'target', 'target tree unconnected', 'distractor', 'distractor unconnected'});end
            plot(tb,smooth(ttree_t_avg{2,ct},smoothness),'b', 'LineWidth',2);
            plot(tb,smooth(ttree_d_avg{2,ct},smoothness),'r', 'LineWidth',2);
            plot(tb,smooth(dtree_t_avg{2,ct},smoothness),'c', 'LineWidth',2);
            plot(tb,smooth(dtree_d_avg{2,ct},smoothness),'m', 'LineWidth',2); xlim([-.4 0.9]);
            if Poss_Cue_Time(ct)>=0
                line([fix_stim_dur+1e-3*Poss_Cue_Time(ct) fix_stim_dur+1e-3*Poss_Cue_Time(ct)], [-10 20])
            end
            ylim([-.2 1.2])
        end
        
        figure(f2)
        filename = ['NeurActiv_withStimDur_Cue_' StimCondStr{stimcond} 'Stim_s' num2str(Sessions(1)) '_' num2str(Sessions(end)) '_chan' num2str(ana_chan(1)) '_' num2str(ana_chan(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        
        %% plot all the cue times for each conditions
        figure(f3);clf
        subplot(2,2,1);hold on; title(['Connected - Target pair'])%
        labels = cell(1,numel(Poss_Cue_Time));
        for ct = 1:numel(Poss_Cue_Time)
            plot(tb,smooth(ttree_t_avg{2,ct},smoothness),'Color', [max(0,.2*ct-.4) max(0,.2*ct-.4) min(1,.4*ct)]);
            labels{ct} = ['Cue at ' num2str(Poss_Cue_Time(ct)) 'ms'];
        end
        xlim([-.2 .9]);
        legend(labels)
        subplot(2,2,2);hold on; title(['Unconnected - Target pair'])%
        for ct = 1:numel(Poss_Cue_Time)
            plot(tb,smooth(ttree_d_avg{2,ct},smoothness),'Color', [min(1,.4*ct) max(0,.2*ct-.4) max(0,.2*ct-.4)]);
        end
        xlim([-.2 .9]);
        subplot(2,2,3);hold on; title(['Connected - Distr pair'])%
        for ct = 1:numel(Poss_Cue_Time)
            plot(tb,smooth(dtree_t_avg{2,ct},smoothness),'Color', [max(0,.2*ct-.4) min(1,.4*ct) min(1,.4*ct)]);
        end
        xlim([-.2 .9]);
        subplot(2,2,4);hold on; title(['Unconnected - Distr pair'])%
        for ct = 1:numel(Poss_Cue_Time)
            plot(tb,smooth(dtree_d_avg{2,ct},smoothness),'Color', [min(1,.4*ct) max(0,.2*ct-.4) min(1,.4*ct)]);
        end
        xlim([-.2 .9]);
        filename = ['NeurActiv_withCond_Cue_' 's' num2str(Sessions(1)) '_' num2str(Sessions(end)) '_chan' num2str(ana_chan(1)) '_' num2str(ana_chan(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        
        %% plot
        %         figure('Position', [100 50 1400 900])
        %         AUC_avg = zeros(1,min(numel(Poss_Cue_Time),6));
        %         for ct = 1:min(numel(Poss_Cue_Time),6)
        %             subplot(2,4,ct);hold on; title(['cue delay:' num2str(Poss_Cue_Time(ct)), ' ms'])
        %             ylim1 = [-2 2.2];
        %             xlabel('time(s)'); ylabel('normalized activity')
        %
        %             ttree_t_avg{2,ct} = ttree_t_avg{2,ct}-mean(ttree_t_avg{2,ct}(tb<baseline(2)));
        %             dtree_t_avg{2,ct} = dtree_t_avg{2,ct}-mean(dtree_t_avg{2,ct}(tb<baseline(2)));
        %             a = plot(tb, smooth(ttree_t_avg{2,ct}, smoothness), 'b'); ylim(ylim1)
        %             [hAx,hLine1,hLine2] = plotyy(tb, smooth(dtree_t_avg{2,ct}, smoothness), tb, smooth(ttree_t_avg{2,ct}-dtree_t_avg{2,ct}, smoothness));
        %             hAx(2).XLim = x_limits;
        %             hAx(1).XLim = x_limits;
        %             hLine1.Color = [0 1 1];
        %             hAx(1).YLim = ylim1;
        %             hAx(2).YLim = [-.6 2];
        %             hLine2.Color=[0 0 0];
        %             axes(hAx(2));line(hAx(2).XLim,[0 0])
        %             if ct==1
        %                 legend([a, hLine1,hLine2], { 'target', 'connected distractor', 'difference'})
        %             end
        %
        %             if Poss_Cue_Time(ct)<0
        %                 p=patch([0 integ_dur integ_dur 0],[-10 -10 20 20],'r');
        %                 t_relev = find(tb>(0) & tb<(integ_dur));
        %             else
        %                 p=patch([fix_stim_dur, fix_stim_dur+integ_dur, fix_stim_dur+integ_dur, fix_stim_dur],[-10 -10 20 20],'r');
        %                 t_relev = find(tb>(fix_stim_dur) & tb<(fix_stim_dur +integ_dur));
        %             end
        %             set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
        %             diff =abs(smooth(ttree_t_avg{2,ct}, smoothness)-smooth(dtree_t_avg{2,ct}, smoothness));
        %             diff = diff(t_relev);diff(diff<0)=0;
        %             AUC_avg(ct) = sum(diff);
        %             if Poss_Cue_Time(ct)>=0
        %                 line([fix_stim_dur+1e-3*Poss_Cue_Time(ct) fix_stim_dur+1e-3*Poss_Cue_Time(ct)], [-10 20])
        %             end
        %             %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
        %         end
        %         %             subplot(2,4,ct+1);hold on; title('no cue')
        %         %             ttree_t_avg{2} = ttree_t_avg{2}-mean(ttree_t_avg{2}(tb<baseline(2)));
        %         %             dtree_t_avg{2} = dtree_t_avg{2}-mean(dtree_t_avg{2}(tb<baseline(2)));
        %         %             b = plot(tb, smooth(ttree_t_avg{2}, smoothness),'Color', [0 0 .7]);  ylim(ylim1);
        %         %             [hAx2,hLine12,hLine22] = plotyy(tb, smooth(dtree_t_avg{2}, smoothness),tb, smooth(ttree_t_avg{2}-dtree_t_avg{2}, smoothness)); xlim([-.2 0.9]);
        %         %             hLine12.Color = [0 .7 .7];
        %         %             hLine22.Color=[0 0 0];
        %         %             hAx2(2).XLim = x_limits;
        %         %             hAx2(2).YLim = [-.6 2];
        %         %             hAx2(1).YLim = ylim1;
        %         %             axes(hAx2(2));line(hAx2(2).XLim,[0 0])
        %         %             p=patch([fix_stim_dur fix_stim_dur+integ_dur fix_stim_dur+integ_dur fix_stim_dur],[-10 -10 20 20],'r');
        %         %             set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
        %         %             legend([b ,hLine12,hLine22], {'target', 'connected distractor', 'difference'})
        %
        %         subplot(2,4,8);hold on;
        %         %          AUC_avg = 0;
        % %         t_relev = find(tb>(fix_stim_dur) & tb<(fix_stim_dur +integ_dur));
        %
        %
        %
        %         hold on; plot(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_avg, 'g')
        %
        %         %     legend({'no mask', 'masked'})
        %         xlabel('cue delay(ms)')
        %         title('Area Under the difference curve from stim end to stim end+300ms')
        %         e = std(AUC(:,good_chann),1,2)'/2;
        %         [hAx,hLine1,hLine2] = plotyy(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_avg,[Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6))],  [perf(2,1:min(numel(Poss_Cue_Time),6))./sum(perf(:,1:min(numel(Poss_Cue_Time),6))) ]);
        %         hold on; errorbar(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_avg,e, 'g')
        %         hLine2.Color = [1 0 0];
        %         hLine1.Color = [0 .8 0];
        %         hAx(1).YLim = [0 40];
        %         hAx(2).YLim = [.5 1];
        %         legend([hLine1;hLine2],'neuronal','behavior','Location','northeast');
        %
        %         filename = ['Diff_ConnectedUnconn_' StimCondStr{stimcond} 'Stim_s' num2str(Sessions(1)) '_' num2str(Sessions(end)) '_chan' num2str(ana_chan(1)) '_' num2str(ana_chan(end))];
        %         set(gcf, 'PaperPositionMode', 'auto')
        %         print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        %
        %
        %% plot the difference between connected and unconnected for target and distr pair
        figure('Position', [100 50 1400 900])
        AUC_avg = zeros(1,min(numel(Poss_Cue_Time),6));
        for ct = 1:min(numel(Poss_Cue_Time),6)
            subplot(2,4,ct);hold on; title(['cue delay:' num2str(Poss_Cue_Time(ct)), ' ms'])
            ylim1 = [-1.3 1];
            xlabel('time(s)'); ylabel('normalized activity')
            ttree_t_avg{2,ct} = ttree_t_avg{2,ct}-mean(ttree_t_avg{2,ct}(tb<baseline(2)));
            ttree_d_avg{2,ct} = ttree_d_avg{2,ct}-mean(ttree_d_avg{2,ct}(tb<baseline(2)));
            dtree_t_avg{2,ct} = dtree_t_avg{2,ct}-mean(dtree_t_avg{2,ct}(tb<baseline(2)));
            dtree_d_avg{2,ct} = dtree_d_avg{2,ct}-mean(dtree_d_avg{2,ct}(tb<baseline(2)));
            a = plot(tb, smooth(ttree_t_avg{2,ct}, smoothness)-smooth(ttree_d_avg{2,ct}, smoothness), 'g'); ylim(ylim1)
            [hAx,hLine1,hLine2] = plotyy(tb, smooth(dtree_t_avg{2,ct}, smoothness)-smooth(dtree_d_avg{2,ct}, smoothness), tb, smooth(ttree_t_avg{2,ct}-ttree_d_avg{2,ct}-dtree_t_avg{2,ct}+dtree_d_avg{2,ct}, smoothness));
            hAx(2).XLim = x_limits;
            hAx(1).XLim = x_limits;
            hLine1.Color = [1 .6 0];
            hAx(1).YLim = ylim1;
            hAx(2).YLim = [-.6 2];
            hLine2.Color=[0 0 0];
            axes(hAx(2));line(hAx(2).XLim,[0 0])
            if ct==1
                legend([a, hLine1,hLine2], { 'target pair', ' distractor pair', 'difference'})
            end
            
            if Poss_Cue_Time(ct)<0
                p=patch([0.1 0.1+integ_dur .1+integ_dur 0.1],[-10 -10 20 20],'r');
                t_relev = find(tb>(0.1) & tb<(0.1+integ_dur));
            else
                p=patch([fix_stim_dur+Poss_Cue_Time(ct)*1e-3, fix_stim_dur+Poss_Cue_Time(ct)*1e-3+integ_dur, fix_stim_dur+Poss_Cue_Time(ct)*1e-3+integ_dur, fix_stim_dur+Poss_Cue_Time(ct)*1e-3],[-10 -10 20 20],'r');
                t_relev = find(tb>(fix_stim_dur+Poss_Cue_Time(ct)*1e-3) & tb<(fix_stim_dur+Poss_Cue_Time(ct)*1e-3 +integ_dur));
            end
            set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
            diff =(smooth(ttree_t_avg{2,ct}-ttree_d_avg{2,ct}, smoothness)-smooth(dtree_t_avg{2,ct}-dtree_d_avg{2,ct}, smoothness));
            diff = diff(t_relev);%diff(diff<0)=0;
            AUC_avg(ct) = sum(diff);
            if Poss_Cue_Time(ct)>=0
                line([fix_stim_dur+1e-3*Poss_Cue_Time(ct) fix_stim_dur+1e-3*Poss_Cue_Time(ct)], [-10 20])
            end
        end
        
        subplot(2,4,8);hold on;
        hold on; plot(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_avg, 'g')
        
        xlabel('cue delay(ms)')
        title('Area Under the difference curve from stim end to fp green')
        e = std(AUC(:,good_chann),1,2)'/2;
        [hAx,hLine1,hLine2] = plotyy(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_avg,[Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)) ],  [perf(2,1:min(numel(Poss_Cue_Time),6))./sum(perf(:,1:min(numel(Poss_Cue_Time),6))) ]);
        hold on; errorbar(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_avg,e, 'g')
        hLine2.Color = [1 0 0];
        hLine1.Color = [0 .8 0];
        hAx(1).YLim = [-10 45];
        hAx(2).YLim = [.5 1];
        legend([hLine1;hLine2],'neuronal','behavior','Location','northeast');
        
        filename = ['Diff_InfoTargNonTarg_' StimCondStr{stimcond} 'Stim_s' num2str(Sessions(1)) '_' num2str(Sessions(end)) '_chan' num2str(ana_chan(1)) '_' num2str(ana_chan(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        
        
        neur_decay = [AUC_avg AUC_avg];
        cue_times = [Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)) NaN];
        behav_dacay = [perf(2,1:min(numel(Poss_Cue_Time),6))./sum(perf(:,1:min(numel(Poss_Cue_Time),6)))];
        filename = ['Monty_varCue_NeurBehav_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
        save(filename, 'neur_decay', 'behav_dacay', 'cue_times')
        
        
        %% all diff curves on same graph
        figure;hold on
        smoothness = 40;
        col = colormap('hot');
        col = col(1:10:end,:);
        for ct = 1:min(numel(Poss_Cue_Time),6)
            plot(tb, smooth(ttree_t_avg{2,ct}-ttree_d_avg{2,ct}-dtree_t_avg{2,ct}+dtree_d_avg{2,ct}, smoothness), 'Color', col(ct,:))
        end
        xlim([-.2 0.9]);
%         plot(tb, smooth(ttree_t_avg{2}-dtree_t_avg{2}, smoothness), 'k', 'LineWidth',2); 
        
    end
end


