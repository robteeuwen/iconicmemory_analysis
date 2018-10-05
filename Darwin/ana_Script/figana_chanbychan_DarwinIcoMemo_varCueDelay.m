% analysis function for the cueing of iconic memory  
% This script can be used for the initial cueing  experiment
% (s73-95) 
% written by Catherine Wacongne 

close all
clear all
dbstop if error
plot_stim_dur =1;
plot_mask = 1;
plot_chann= 1;

fix_stim_dur = 0.06;
integ_dur = 0.6;
smoothness = 30;

[extractdir, sessions, rawdir, SaveDir] = infoDir_DarwinIcoMemo;
info = Log_DarwinIcoMemo;


load([extractdir, sessions{end} '\EVT_', info(end).Tankname,'_block-',num2str(info(end).goodblocs(1))])


SF = EVENT.strms(1).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);


Cond = cell(96,2,2,2);
ci = zeros(96,2);

baseline = [-.5 -.2];
x_limits = [-.5 1.2];
%73 to 76: precue with disapearing segment in RF 77:78
%83-86: session with dark distr
%87-90: session with yellow distr
%94 - 95 sessions with pre-post cue
%96-101: new sessions after retraining - recorded by Rob, based on Monty's paradigm
Sessions = 97:101;%
aregood = zeros(max(Sessions),48);
for n = 25:48
    
    % Extract the trials from all blocs and conditions for each channel
    
    AllTrials = [];
    AllCond = [];
    days =[];
    for s = Sessions%
        
        for bloc = 1:numel(info(s).goodblocs)
            %Read in the extracted data form the mapped channel
            load([extractdir sessions{s},'\Xtract_',info(s).Tankname,'_Block-',num2str(info(s).goodblocs(bloc)),'_',num2str(n)]);
            e = Env{1};
            load([rawdir,'ICOmemo_' info(s).Tankname '_B',num2str(info(s).goodblocs(bloc))])
            trials =find(LOG.target_presented>0);
            equal = (size(e,2) == length(trials));
            if size(e,2)<length(trials)
                trials = trials(1:size(e,2));
            elseif size(e,2)>length(trials)
                disp(['skipping bloc ' num2str(bloc) ' in session ' num2str(s)])
                continue
            end
            %             if ~equal
            %                 disp(['skipping bloc ' num2str(bloc) ' in session ' num2str(s)])
            %                 continue
            %             end
            isdrawn4_b = zeros(1,numel(trials));
            connected_b = zeros(1,numel(trials));
            for i =1: numel(trials)
                a = find(LOG.drawn_skeleton{trials(i)}==3);
                if ~isempty(a)
                    isdrawn4_b(i) = a;
                    connected_b(i) = LOG.drawn_conn{trials(i)}(a);
                end
                
            end
            days = [days, s*ones(1,size(e,2))];
            if isempty(AllTrials)
                AllTrials=e;
                isdrawn4 = isdrawn4_b;
                connected = connected_b;
                targ = LOG.targ_num(trials);
                Hit = LOG.correct(trials)*2 + LOG.congruent_error(trials)+LOG.incongruent_error(trials) + LOG.incongruent_hit(trials);
                %Hit = LOG.correct(trials)+1;
                stim_dur = LOG.Stim_dur(trials);
                cue_time = LOG.CueTime(trials);
                Poss_Cue_Time = unique(cue_time);
                mask = LOG.ismask(trials)+1;
                mask = mask.*(LOG.CueDuration(trials)==100) +1;
                
                %                 AllCond = MAT;
            else
                AllTrials = [AllTrials e];
                %                 AllCond = [AllCond;MAT];
                isdrawn4 = [isdrawn4 isdrawn4_b];
                connected = [connected connected_b];
                targ = [targ LOG.targ_num(trials)];
                Hit = [Hit (LOG.correct(trials)*2  + LOG.congruent_error(trials)+LOG.incongruent_error(trials) + LOG.incongruent_hit(trials))];
                
                %                 Hit = [Hit LOG.correct(trials)+1];
                stim_dur = [stim_dur LOG.Stim_dur(trials)];
                cue_time = [cue_time LOG.CueTime(trials)];
                Poss_Cue_Time = [Poss_Cue_Time unique(cue_time)];
                mask_d = LOG.ismask(trials)+1;
                mask = [mask mask_d.*(LOG.CueDuration(trials)==100)+1];
            end
            
        end
    end
    Poss_Cue_Time = unique(Poss_Cue_Time);
    % do some cleaning
    %     ev = reshape(AllTrials,size(AllTrials,1)*size(AllTrials,2),1);
    %     ev(abs(zscore(ev))>5) = NaN;
    %     AllTrials = reshape(ev,size(AllTrials,1),size(AllTrials,2)); figure;subplot(1,2,1);imagesc(AllTrials);
    %
    %     %Z-score the trials
    %     T = nanmean(AllTrials);
    %     f1 = find(abs(zscore(T))>3);
    %     f2 = find(max(abs(AllTrials))>0.1);
    %     AllTrials(:,[f1,f2]) = NaN(size(AllTrials,1),length(f1)+length(f2));subplot(1,2,2);imagesc(AllTrials);pause
    %     keyboard
    bf = find(tb<baseline(2));
    pf = find(tb > 0.03 & tb < 0.1);
    pf2 = find(tb > 0.05 & tb < 0.09);
    hitNorm =(targ==5 & Hit==2);
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
    for h = 1:2
        if plot_stim_dur
            
            % average per condition
            ind = find(isdrawn4==0 & Hit==h  & mask ==1 & stim_dur<100);%& stim_dur == Poss_Cue_Time(sd) & mask == m
            background{n,h}  = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(targ==5 & Hit==h  & mask ==1 & stim_dur<100);
            ttree_t{n,h} = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(targ==6 & Hit==h  & mask ==1 & stim_dur<100);
            ttree_d{n,h} = nanmean(AllTrials(:,ind),2);
            
            ind = find(connected==1 & targ~=5 & Hit==h & mask ==1 & stim_dur<100);
            dtree_t{n,h} = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(connected==2 & targ~=6 & Hit==h & mask ==1 & stim_dur<100);
            dtree_d{n,h} = nanmean(AllTrials(:,ind),2);
            
            ind = find( Hit==h  & mask ==1 & stim_dur<100 );
            perf(h) = numel(ind);
            
            
            ind1 = find(connected==1  & mask ==1 & Hit==2 & stim_dur<100);%
            t_connected{n} =  nanmean(AllTrials(:,ind1),2);
            ind2 = find(connected==2  & mask ==1 & Hit==2 & stim_dur<100);%
            t_unconnected{n} = nanmean(AllTrials(:,ind2),2);
            %             figure; hist(nanmean(AllTrials(pf2,ind1)),20); hold on; hist(nanmean(AllTrials(pf2,ind2)),20);title(num2str(n));pause(0.1);
            [hyp, ismodul(n),ci(n,:),stats] = ttest2(nanmean(AllTrials(pf2,ind1)),nanmean(AllTrials(pf2,ind2))) ;
            if plot_mask
                for sd = 1:numel(Poss_Cue_Time)
                    ind = find(isdrawn4==0 & Hit==h & cue_time == Poss_Cue_Time(sd) & mask ==3 & stim_dur<100);%& stim_dur == Poss_Cue_Time(sd) & mask == m
                    background_m{n,h,sd}  = nanmean(AllTrials(:,ind),2);
                    
                    
                    ind = find(targ==5 & Hit==h & cue_time == Poss_Cue_Time(sd) & mask ==3 & stim_dur<100);
                    ttree_t_m{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    
                    
                    ind = find(targ==6 & Hit==h & cue_time == Poss_Cue_Time(sd) & mask ==3 & stim_dur<100);
                    ttree_d_m{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    ind = find(connected==1 & targ~=5 & Hit==h & cue_time == Poss_Cue_Time(sd)& mask ==3  & stim_dur<100);
                    dtree_t_m{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    
                    ind = find(connected==2 & targ~=6 & Hit==h & cue_time == Poss_Cue_Time(sd)& mask ==3 & stim_dur<100 );
                    dtree_d_m{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    ind = find( Hit==h & cue_time == Poss_Cue_Time(sd) & mask ==3 & stim_dur<100);
                    perf_m(h,sd) = numel(ind);
                    
                    
                    ind = find(connected==1 & cue_time == Poss_Cue_Time(sd) & mask ==3 & Hit==2 & stim_dur<100);%
                    t_connected_m{n,sd} =  nanmean(AllTrials(:,ind),2);
                    ind = find(connected==2 & cue_time == Poss_Cue_Time(sd) & mask ==3 & Hit==2 & stim_dur<100);%
                    t_unconnected_m{n,sd} = nanmean(AllTrials(:,ind),2);
                    
                end
                
            end
        else
            ind = find(isdrawn4==0 & Hit==h & mask ==1  & stim_dur<100);%& stim_dur == Poss_Cue_Time(sd) & mask == m
            background{n,h}  = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(targ==5 & Hit==h & mask ==1  & stim_dur<100);
            ttree_t{n,h} = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(targ==6 & Hit==h & mask ==1  & stim_dur<100);
            ttree_d{n,h} = nanmean(AllTrials(:,ind),2);
            
            ind = find(connected==1 & targ~=5 & Hit==h & mask ==1 & stim_dur<100);
            dtree_t{n,h} = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(connected==2 & targ~=6 & Hit==h & mask ==1 & stim_dur<100);
            dtree_d{n,h} = nanmean(AllTrials(:,ind),2);
            
            if plot_mask
                ind = find(isdrawn4==0 & Hit==h & mask ==3  & stim_dur<100);%& stim_dur == Poss_Cue_Time(sd) & mask == m
                background_m{n,h}  = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(targ==5 & Hit==h & mask ==3  & stim_dur<100);
                ttree_t_m{n,h} = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(targ==6 & Hit==h & mask ==3  & stim_dur<100);
                ttree_d_m{n,h} = nanmean(AllTrials(:,ind),2);
                
                ind = find(connected==1 & targ~=5 & Hit==h & mask ==3 & stim_dur<100);
                dtree_t_m{n,h} = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(connected==2 & targ~=6 & Hit==h & mask ==3 & stim_dur<100);
                dtree_d_m{n,h} = nanmean(AllTrials(:,ind),2);
                
            end
            
        end
    end
    
end

%% plot the conditions

close all
if plot_chann
    %     f1 =figure('Position', [100 50 1400 900]);
    f2 =figure('Position', [100 50 1400 900]);
    %n,h,sd,m
    AUC = zeros(1,48);
    AUC_m = zeros(min(numel(Poss_Cue_Time),6),48);
    for n = 25:48
        
        
        if plot_stim_dur
            
            %             figure(f1);clf
            figure(f2);clf
            for sd = 1:min(numel(Poss_Cue_Time),6)
                %                  figure(f1);subplot(2,4,sd);hold on
                %                 plot(tb,smooth(background{n,1,sd},smoothness),'k');
                %                 plot(tb,smooth(ttree_t{n,1,sd},smoothness),'b');
                %                 %         plot(tb,smooth(ttree_d{n,1,sd},smoothness),'r');
                %                 %         plot(tb,smooth(dtree_t{n,1,sd},smoothness),'c');
                %                 %         plot(tb,smooth(dtree_d{n,1,sd},smoothness),'m');
                %                 %         legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'})
                %                 plot(tb,smooth(background{n,2,sd},smoothness),'k', 'LineWidth',2);
                %                 plot(tb,smooth(ttree_t{n,2,sd},smoothness),'b', 'LineWidth',2); xlim([-.2 1.1]);
                %                 %         plot(tb,smooth(ttree_d{n,2,sd},smoothness),'r', 'LineWidth',2);
                %                 plot(tb,smooth(dtree_t{n,2,sd},smoothness),'c', 'LineWidth',2);
                %                 %         plot(tb,smooth(dtree_d{n,2,sd},smoothness),'m', 'LineWidth',2);
                %                 title(num2str(n))
                %
                
                
                figure(f2)
                subplot(2,4,sd);hold on; title([num2str(Poss_Cue_Time(sd)), ' ms, chan ' num2str(n)])
                
                if plot_mask
                    plot(tb, smooth(ttree_t_m{n,2,sd}, smoothness), 'b', 'LineWidth',1); try  ylim([min(smooth(ttree_t_m{n,2,sd},smoothness))-(max(smooth(ttree_t_m{n,2,sd},smoothness))-min(smooth(ttree_t_m{n,2,sd},smoothness))),  max(smooth(ttree_t_m{n,2,sd},smoothness))+.2*(max(smooth(ttree_t_m{n,2,sd},smoothness))-min(smooth(ttree_t_m{n,2,sd},smoothness)))]);end
                    [hAx,hLine1,hLine2] = plotyy(tb, smooth(dtree_t_m{n,2,sd}, smoothness), tb, smooth(ttree_t_m{n,2,sd}-dtree_t_m{n,2,sd}, smoothness));
                    diffcu =smooth(ttree_t_m{n,2,sd}-dtree_t_m{n,2,sd},smoothness);
                    hAx(2).XLim = x_limits;
                    hAx(1).XLim = x_limits;
                    hLine1.LineWidth = 1;
                    hLine1.Color = [1 0 0];
                    try hAx(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
                    hLine2.Color=[0 1 0];
                    hLine2.LineWidth = 1;
                    
                    %             plot(tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, smoothness), 'g', 'LineWidth',2);
                end
                %                 plot(tb, smooth(t_connected{n,sd}, smoothness),'b')
                %                 [hAx2,hLine12,hLine22] = plotyy(tb, smooth(t_unconnected{n,sd}, smoothness),tb, smooth(t_connected{n,sd}-t_unconnected{n,sd}, smoothness)); xlim([-.2 0.9]);
                %                 hLine12.Color = [1 0 0];
                %                 hLine22.Color=[0 0 0];
                %                 hAx2(2).XLim = [-.2 0.9];
                axes(hAx(2));line(hAx(2).XLim,[0 0]);
                p=patch([fix_stim_dur fix_stim_dur+integ_dur fix_stim_dur+integ_dur fix_stim_dur],[-10 -10 20 20],'r');
                set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
                %                 try;hAx(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
                %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
                
                
                
                
            end
            
            subplot(2,4,sd+1);hold on; title('no cue')
            plot(tb, smooth(ttree_t{n,2}, smoothness),'b');  try  ylim([min(smooth(ttree_t{n,2},smoothness))-(max(smooth(ttree_t{n,2},smoothness))-min(smooth(ttree_t{n,2},smoothness))),  max(smooth(ttree_t{n,2},smoothness))+.2*(max(smooth(ttree_t{n,2},smoothness))-min(smooth(ttree_t{n,2},smoothness)))]);end
            [hAx,hLine1,hLine2] = plotyy(tb, smooth(dtree_t{n,2}, smoothness),tb, smooth(ttree_t{n,2}-dtree_t{n,2}, smoothness)); xlim([-.2 0.9]);
            hLine1.Color = [1 0 0];
            hLine2.Color=[0 0 0];
            hAx(1).XLim = x_limits;
            hAx(2).XLim = x_limits;
            axes(hAx(2));line(hAx(2).XLim,[0 0]);
            p=patch([fix_stim_dur fix_stim_dur+integ_dur fix_stim_dur+integ_dur fix_stim_dur],[-10 -10 20 20],'r');
            set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
            try hAx(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
            
            %             AUC = 0;
            t_relev = find(tb>(fix_stim_dur) & tb<(fix_stim_dur +integ_dur));
            diff = smooth(ttree_t{n,2},smoothness)- smooth(dtree_t{n,2},smoothness);
            diff = diff(t_relev);diff(diff<0)=0;
            %             i = find(diff<0, 1);diff(i:end)=0;
            AUC(n) = sum(diff);
            
            if plot_mask
                %                 AUC_m = zeros(1,6);
                for sd = 1:min(numel(Poss_Cue_Time),6)
                    
                    diff_m =smooth(ttree_t_m{n,2,sd}, smoothness)-smooth(dtree_t_m{n,2,sd}, smoothness);
                    diff_m = diff_m(t_relev);diff_m(diff_m<0)=0;
                    %                     i = find(diff_m<0, 1);diff_m(i:end)=0;
                    AUC_m(sd,n) = sum(diff_m);
                end
                subplot(2,4,8);hold on; plot(400, AUC(n), 'bo'); plot(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_m(:,n), 'r')
            end
            
            legend({'no cue', 'cued'})
            xlabel('cue delay')
            pause(0.2)
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-dpng',[SaveDir 'ConnUnconn_Chan' num2str(n) '_s' num2str(days(1)) '_' num2str(days(end)) '.png'])
            
        else
            clf
            
            hold on
            plot(tb,smooth(background{n,1},smoothness),'k');
            plot(tb,smooth(ttree_t{n,1},smoothness),'b');
            plot(tb,smooth(ttree_d{n,1},smoothness),'r');
            plot(tb,smooth(dtree_t{n,1},smoothness),'c');
            plot(tb,smooth(dtree_d{n,1},smoothness),'m');
            legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'})
            plot(tb,smooth(background{n,2},smoothness),'k', 'LineWidth',2);
            plot(tb,smooth(ttree_t{n,2},smoothness),'b', 'LineWidth',2);
            plot(tb,smooth(ttree_d{n,2},smoothness),'r', 'LineWidth',2);
            plot(tb,smooth(dtree_t{n,2},smoothness),'c', 'LineWidth',2);
            plot(tb,smooth(dtree_d{n,2},smoothness),'m', 'LineWidth',2); xlim([-.2 1.1]);
            title(num2str(n)); pause
            
            
            
        end
        
        
    end
    
    
    
    
end
%%
% cand = find(ismodul<0.001);
% cand = find(ci(:,1)>0.1)';
% good_chann = cand(cand>72);%[75:76 78 79 81 83:87 89 90 92:94 96];%74:96;%[75 81 82 83 85:89 91:96];%[74 75 78 83 86 90:96];
good_chann = find(SNR>.5);
if plot_stim_dur
    figure; plot(400,  perf(2)./sum(perf), 'bo')
    fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
    title('Behavioral performance on the iconic memory task');
    xlabel('Stimulus-cue interval (ms)')
    if plot_mask
        hold on;  plot(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)),  perf_m(2,1:min(numel(Poss_Cue_Time),6))./sum(perf_m(:,1:min(numel(Poss_Cue_Time),6))),'r')
        legend('no cue', 'cued')
        
    end
    
    filename = 'performances';
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    %%
    f1 = figure('Position', [100 50 1400 900]);
    if plot_mask; f2 = figure('Position', [100 50 1400 900]);end
    for sd = 1:min(numel(Poss_Cue_Time),6)
        t_connected_avg = 0*background{25,1};n_t_connected_avg=0;
        t_unconnected_avg = 0*background{25,1};n_t_unconnected_avg=0;
        for h = 1:2
            background_avg{h} = 0*background{25,1};n_background_avg(h)=0;
            ttree_t_avg{h} = 0*background{25,1};n_ttree_t_avg(h)=0;
            ttree_d_avg{h} = 0*background{25,1};n_ttree_d_avg(h)=0;
            dtree_t_avg{h} = 0*background{25,1};n_dtree_t_avg(h)=0;
            dtree_d_avg{h} = 0*background{25,1};n_dtree_d_avg(h)=0;
        end
        if plot_mask
            t_connected_m_avg{sd} = 0*background{25,1};n_t_connected_m_avg(sd)=0;
            t_unconnected_m_avg{sd} = 0*background{25,1};n_t_unconnected_m_avg(sd)=0;
            for h = 1:2
                background_m_avg{h,sd} = 0*background{25,1};n_background_m_avg(h,sd)=0;
                ttree_t_m_avg{h,sd} = 0*background{25,1};n_ttree_t_m_avg(h,sd)=0;
                ttree_d_m_avg{h,sd} = 0*background{25,1};n_ttree_d_m_avg(h,sd)=0;
                dtree_t_m_avg{h,sd} = 0*background{25,1};n_dtree_t_m_avg(h,sd)=0;
                dtree_d_m_avg{h,sd} = 0*background{25,1};n_dtree_d_m_avg(h,sd)=0;
            end
        end
        for n = good_chann
            
            try t_connected_avg = nansum([t_connected_avg ,t_connected{n}],2);end; if ~isnan(sum(t_connected{n}));n_t_connected_avg = n_t_connected_avg+1;end
            try t_unconnected_avg =  nansum([t_unconnected_avg,t_unconnected{n}],2);end;if ~isnan(sum(t_unconnected{n}));n_t_unconnected_avg = n_t_unconnected_avg+1;end
            for h = 1:2
                try background_avg{h} =  nansum([background_avg{h} ,background{n,h}],2);end;if ~isnan(sum(background{n,h}));n_background_avg(h) = n_background_avg(h)+1;end
                try ttree_t_avg{h} =  nansum([ttree_t_avg{h} ,ttree_t{n,h}],2);end ;if ~isnan(sum(ttree_t{n,h}));n_ttree_t_avg(h) = n_ttree_t_avg(h)+1;end
                try ttree_d_avg{h} =  nansum([ttree_d_avg{h} ,ttree_d{n,h} ],2);end;if ~isnan(sum(ttree_d{n,h}));n_ttree_d_avg(h) = n_ttree_d_avg(h)+1;end
                
                try dtree_t_avg{h} =  nansum([dtree_t_avg{h},dtree_t{n,h}],2);end;if ~isnan(sum(dtree_t{n,h}));n_dtree_t_avg(h) = n_dtree_t_avg(h)+1;end
                try dtree_d_avg{h} =  nansum([dtree_d_avg{h},dtree_d{n,h}],2);end;if ~isnan(sum(dtree_d{n,h}));n_dtree_d_avg(h) = n_dtree_d_avg(h)+1;end
            end
            if plot_mask
                try t_connected_m_avg{sd} =  nansum([t_connected_m_avg{sd} ,t_connected_m{n,sd}],2);end;if ~isnan(sum(t_connected_m{n,sd}));n_t_connected_m_avg(sd) = n_t_connected_m_avg(sd)+1;end
                try t_unconnected_m_avg{sd} =  nansum([t_unconnected_m_avg{sd} ,t_unconnected_m{n,sd}],2);end;if ~isnan(sum(t_unconnected_m{n,sd}));n_t_unconnected_m_avg(sd) = n_t_unconnected_m_avg(sd)+1;end
                for h = 1:2
                    try background_m_avg{h,sd} =  nansum([background_m_avg{h,sd} ,background_m{n,h,sd}],2);end;if ~isnan(sum(background_m{n,h,sd}));n_background_m_avg(h,sd) = n_background_m_avg(h,sd)+1;end
                    try  ttree_t_m_avg{h,sd} =  nansum([ttree_t_m_avg{h,sd} ,ttree_t_m{n,h,sd}],2);end;if ~isnan(sum(ttree_t_m{n,h,sd}));n_ttree_t_m_avg(h,sd) = n_ttree_t_m_avg(h,sd)+1;end
                    try ttree_d_m_avg{h,sd} =  nansum([ttree_d_m_avg{h,sd} ,ttree_d_m{n,h,sd} ],2);end;if ~isnan(sum(ttree_d_m{n,h,sd}));n_ttree_d_m_avg(h,sd) = n_ttree_d_m_avg(h,sd)+1;end
                    
                    try dtree_t_m_avg{h,sd} =  nansum([dtree_t_m_avg{h,sd},dtree_t_m{n,h,sd}],2);end;if ~isnan(sum(dtree_t_m{n,h,sd}));n_dtree_t_m_avg(h,sd) = n_dtree_t_m_avg(h,sd)+1;end
                    try dtree_d_m_avg{h,sd} =  nansum([dtree_d_m_avg{h,sd},dtree_d_m{n,h,sd}],2);end;if ~isnan(sum(dtree_d_m{n,h,sd}));n_dtree_d_m_avg(h,sd) = n_dtree_d_m_avg(h,sd)+1;end
                end
            end
        end
        
        t_connected_avg = t_connected_avg/n_t_connected_avg;
        t_unconnected_avg = t_unconnected_avg/n_t_unconnected_avg;
        for h = 1:2
            background_avg{h} = background_avg{h}/n_background_avg(h);
            ttree_t_avg{h} = ttree_t_avg{h}/n_ttree_t_avg(h);
            ttree_d_avg{h} = ttree_d_avg{h}/n_ttree_d_avg(h);
            dtree_t_avg{h} = dtree_t_avg{h}/n_dtree_t_avg(h);
            dtree_d_avg{h} = dtree_d_avg{h}/n_dtree_d_avg(h);
        end
        if plot_mask
            t_connected_m_avg{sd} = t_connected_m_avg{sd}/n_t_connected_m_avg(sd);
            t_unconnected_m_avg{sd} = t_unconnected_m_avg{sd}/n_t_unconnected_m_avg(sd);
            for h = 1:2
                background_m_avg{h,sd} = background_m_avg{h,sd}/n_background_m_avg(h,sd);
                ttree_t_m_avg{h,sd} = ttree_t_m_avg{h,sd}/n_ttree_t_m_avg(h,sd);
                ttree_d_m_avg{h,sd} = ttree_d_m_avg{h,sd}/n_ttree_d_m_avg(h,sd);
                dtree_t_m_avg{h,sd} = dtree_t_m_avg{h,sd}/n_dtree_t_m_avg(h,sd);
                dtree_d_m_avg{h,sd} = dtree_d_m_avg{h,sd}/n_dtree_d_m_avg(h,sd);
            end
        end
        
        
        figure(f1);title('no cue');hold on; %%subplot(2,4,sd)
%         plot(tb,smooth(background_avg{1},smoothness),'k');
%         plot(tb,smooth(ttree_t_avg{1},smoothness),'b');
%         plot(tb,smooth(ttree_d_avg{1},smoothness),'r');
%         plot(tb,smooth(dtree_t_avg{1},smoothness),'c');
%         plot(tb,smooth(dtree_d_avg{1},smoothness),'m');
        if sd == 1; legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'});end
        plot(tb,smooth(background_avg{2},smoothness),'k', 'LineWidth',2);
        plot(tb,smooth(ttree_t_avg{2},smoothness),'b', 'LineWidth',2);
        plot(tb,smooth(ttree_d_avg{2},smoothness),'r', 'LineWidth',2);
        plot(tb,smooth(dtree_t_avg{2},smoothness),'c', 'LineWidth',2);
        plot(tb,smooth(dtree_d_avg{2},smoothness),'m', 'LineWidth',2); xlim([-.3 1.2]);
        if plot_mask
            figure(f2);subplot(2,4,sd);hold on; title([num2str(Poss_Cue_Time(sd)), ' ms, cued at 20ms'])%
%             plot(tb,smooth(background_m_avg{1,sd},smoothness),'k');
%             plot(tb,smooth(ttree_t_m_avg{1,sd},smoothness),'b');
%             plot(tb,smooth(ttree_d_m_avg{1,sd},smoothness),'r');
%             plot(tb,smooth(dtree_t_m_avg{1,sd},smoothness),'c');
%             plot(tb,smooth(dtree_d_m_avg{1,sd},smoothness),'m');
            if sd == 1; legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'});end
            plot(tb,smooth(background_m_avg{2,sd},smoothness),'k', 'LineWidth',2);
            plot(tb,smooth(ttree_t_m_avg{2,sd},smoothness),'b', 'LineWidth',2);
            plot(tb,smooth(ttree_d_m_avg{2,sd},smoothness),'r', 'LineWidth',2);
            plot(tb,smooth(dtree_t_m_avg{2,sd},smoothness),'c', 'LineWidth',2);
            plot(tb,smooth(dtree_d_m_avg{2,sd},smoothness),'m', 'LineWidth',2); xlim([-.3 1.2]);
            
            
        end
    end
    figure(f1)
    filename = 'NeurActiv_withStimDur_noCue';
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    if plot_mask
        figure(f2)
        filename = 'NeurActiv_withStimDur_Cue';
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    end
    %%
    figure('Position', [100 50 1400 900])
    for sd = 1:min(numel(Poss_Cue_Time),6)
        subplot(2,4,sd);hold on; title(['cue delay:' num2str(Poss_Cue_Time(sd)), ' ms'])
        ylim1 = [-2 1.5];
        xlabel('time(s)'); ylabel('normalized activity')
        if plot_mask
            %             t_connected_m_avg{sd} = t_connected_m_avg{sd}-mean(t_connected_m_avg{sd}(tb<0));
            ttree_t_m_avg{2,sd} = ttree_t_m_avg{2,sd}-mean(ttree_t_m_avg{2,sd}(tb<0));
            %             t_unconnected_m_avg{sd} = t_unconnected_m_avg{sd}-mean(t_unconnected_m_avg{sd}(tb<0));
            dtree_t_m_avg{2,sd} = dtree_t_m_avg{2,sd}-mean(dtree_t_m_avg{2,sd}(tb<0));
            a = plot(tb, smooth(ttree_t_m_avg{2,sd}, smoothness), 'b'); ylim(ylim1)
            [hAx,hLine1,hLine2] = plotyy(tb, smooth(dtree_t_m_avg{2,sd}, smoothness), tb, smooth(ttree_t_m_avg{2,sd}-dtree_t_m_avg{2,sd}, smoothness));
            hAx(2).XLim = x_limits;
            hAx(1).XLim = x_limits;
            %             hLine1.LineWidth = 1;
            hLine1.Color = [1 0 0];
            hAx(1).YLim = ylim1;
            hAx(2).YLim = [-.6 2];
            hLine2.Color=[0 1 0];
            %             hLine2.LineWidth = 1;
            axes(hAx(2));line(hAx(2).XLim,[0 0])
            if sd==1
                legend([a, hLine1,hLine2], { 'target', 'connected distractor', 'difference'})
            end
            %             plot(tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, smoothness), 'g', 'LineWidth',2);
        end
        
        p=patch([fix_stim_dur fix_stim_dur+integ_dur fix_stim_dur+integ_dur fix_stim_dur],[-10 -10 20 20],'r');
        set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
        %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
    end
    subplot(2,4,sd+1);hold on; title('no cue')
    ttree_t_avg{2} = ttree_t_avg{2}-mean(ttree_t_avg{2}(tb<0));
    dtree_t_avg{2} = dtree_t_avg{2}-mean(dtree_t_avg{2}(tb<0));
    b = plot(tb, smooth(ttree_t_avg{2}, smoothness),'Color', [0 0 .7]);  ylim(ylim1);
    [hAx2,hLine12,hLine22] = plotyy(tb, smooth(dtree_t_avg{2}, smoothness),tb, smooth(ttree_t_avg{2}-dtree_t_avg{2}, smoothness)); xlim([-.2 0.9]);
    hLine12.Color = [.7 0 0];
    hLine22.Color=[0 0 0];
    hAx2(2).XLim = x_limits;
    hAx2(2).YLim = [-.6 2];
    hAx2(1).YLim = ylim1;
    axes(hAx2(2));line(hAx2(2).XLim,[0 0])
    p=patch([fix_stim_dur fix_stim_dur+integ_dur fix_stim_dur+integ_dur fix_stim_dur],[-10 -10 20 20],'r');
    set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
    legend([b ,hLine12,hLine22], {'target', 'connected distractor', 'difference'})
    
    subplot(2,4,8);hold on;
    %          AUC_avg = 0;
    t_relev = find(tb>(fix_stim_dur) & tb<(fix_stim_dur +integ_dur));
    diff = smooth(ttree_t_avg{2},smoothness)- smooth(dtree_t_avg{2},smoothness);
    diff = diff(t_relev);diff(diff<0)=0;
    %             i = find(diff<0, 1);diff(i:end)=0;
    AUC_avg = sum(diff);
    
    
    
    h1 = errorbar(200, AUC_avg,std(AUC(good_chann)), 'bo');
    if plot_mask
        AUC_m_avg = zeros(1,min(numel(Poss_Cue_Time),6));
        for sd = 1:min(numel(Poss_Cue_Time),6)
            diff_m =smooth(ttree_t_m_avg{2,sd}, smoothness)-smooth(dtree_t_m_avg{2,sd}, smoothness);
            diff_m = diff_m(t_relev);diff_m(diff_m<0)=0;
            %                     i = find(diff_m<0, 1);diff_m(i:end)=0;
            AUC_m_avg(sd) = sum(diff_m);
            %             AUC_m(sd) = sum(smooth(t_connected_m_avg{sd}(t_relev)-t_unconnected_m_avg{sd}(t_relev), smoothness));
        end
        hold on; plot(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_m_avg, 'g')
    end
    %     legend({'no mask', 'masked'})
    xlabel('cue delay(ms)')
    title('Area Under the difference curve from stim end to stim end+300ms')
    e = std(AUC_m(:,good_chann),1,2)'/2;
    [hAx,hLine1,hLine2] = plotyy(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_m_avg,[Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)) 200],  [perf_m(2,1:min(numel(Poss_Cue_Time),6))./sum(perf_m(:,1:min(numel(Poss_Cue_Time),6))) perf(2)/sum(perf)]);
    hold on; errorbar(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_m_avg,e, 'g')
    hLine2.Color = [1 0 0];
    hLine1.Color = [0 .8 0];
    hAx(1).YLim = [0 80];
    hAx(2).YLim = [.5 .9];
    legend([h1;hLine1;hLine2],'no cue','cue','behavior','Location','southeast');
    
    filename = ['Diff_ConnectedUnconn_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    
     %%
            figure('Position', [100 50 1400 900])
            for sd = 1:min(numel(Poss_Cue_Time),6)
                subplot(2,4,sd);hold on; title(['cue delay:' num2str(Poss_Cue_Time(sd)), ' ms'])
                ylim1 = [-2 1.5];
                xlabel('time(s)'); ylabel('normalized activity')
                if plot_mask
                    %             t_connected_m_avg{sd} = t_connected_m_avg{sd}-mean(t_connected_m_avg{sd}(tb<0));
                    ttree_t_m_avg{2,sd} = ttree_t_m_avg{2,sd}-mean(ttree_t_m_avg{2,sd}(tb<0));
                    ttree_d_m_avg{2,sd} = ttree_d_m_avg{2,sd}-mean(ttree_d_m_avg{2,sd}(tb<0));
                    %             t_unconnected_m_avg{sd} = t_unconnected_m_avg{sd}-mean(t_unconnected_m_avg{sd}(tb<0));
                    dtree_t_m_avg{2,sd} = dtree_t_m_avg{2,sd}-mean(dtree_t_m_avg{2,sd}(tb<0));
                    dtree_d_m_avg{2,sd} = dtree_d_m_avg{2,sd}-mean(dtree_d_m_avg{2,sd}(tb<0));
                    a = plot(tb, smooth(ttree_t_m_avg{2,sd}, smoothness)-smooth(ttree_d_m_avg{2,sd}, smoothness), 'b'); ylim(ylim1)
                    [hAx,hLine1,hLine2] = plotyy(tb, smooth(dtree_t_m_avg{2,sd}, smoothness)-smooth(dtree_d_m_avg{2,sd}, smoothness), tb, smooth(ttree_t_m_avg{2,sd}-ttree_d_m_avg{2,sd}-dtree_t_m_avg{2,sd}+dtree_d_m_avg{2,sd}, smoothness));
                    hAx(2).XLim = x_limits;
                    hAx(1).XLim = x_limits;
                    %             hLine1.LineWidth = 1;
                    hLine1.Color = [1 0 0];
                    hAx(1).YLim = ylim1;
                    hAx(2).YLim = [-.6 2];
                    hLine2.Color=[0 1 0];
                    %             hLine2.LineWidth = 1;
                    axes(hAx(2));line(hAx(2).XLim,[0 0])
                    if sd==1
                        legend([a, hLine1,hLine2], { 'target pair', ' distractor pair', 'difference'})
                    end
                    %             plot(tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, smoothness), 'g', 'LineWidth',2);
                end
                
                p=patch([fix_stim_dur fix_stim_dur+integ_dur fix_stim_dur+integ_dur fix_stim_dur],[-10 -10 20 20],'r');
                set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
                %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
            end
            subplot(2,4,sd+1);hold on; title('no cue')
            ttree_t_avg{2} = ttree_t_avg{2}-mean(ttree_t_avg{2}(tb<0));
            dtree_t_avg{2} = dtree_t_avg{2}-mean(dtree_t_avg{2}(tb<0));
            b = plot(tb, smooth(ttree_t_avg{2}, smoothness),'Color', [0 0 .7]);  ylim(ylim1);
            [hAx2,hLine12,hLine22] = plotyy(tb, smooth(dtree_t_avg{2}, smoothness),tb, smooth(ttree_t_avg{2}-dtree_t_avg{2}, smoothness)); xlim([-.2 0.9]);
            hLine12.Color = [.7 0 0];
            hLine22.Color=[0 0 0];
            hAx2(2).XLim = x_limits;
            hAx2(2).YLim = [-.6 2];
            hAx2(1).YLim = ylim1;
            axes(hAx2(2));line(hAx2(2).XLim,[0 0])
            p=patch([fix_stim_dur fix_stim_dur+integ_dur fix_stim_dur+integ_dur fix_stim_dur],[-10 -10 20 20],'r');
            set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
            legend([b ,hLine12,hLine22], {'target pair', 'distractor pair', 'difference'})
            
            subplot(2,4,8);hold on;
            %          AUC_avg = 0;
            t_relev = find(tb>(fix_stim_dur) & tb<(fix_stim_dur +integ_dur));
            diff = abs(smooth(ttree_t_avg{2},smoothness)- smooth(dtree_t_avg{2},smoothness));
            diff = diff(t_relev);diff(diff<0)=0;
            %             i = find(diff<0, 1);diff(i:end)=0;
            AUC_avg = sum(diff);
            
            
            
            h1 = errorbar(200, AUC_avg,std(AUC(good_chann)), 'bo');
            if plot_mask
                AUC_m_avg = zeros(1,min(numel(Poss_Cue_Time),6));
                for sd = 1:min(numel(Poss_Cue_Time),6)
                    diff_m =abs(smooth(ttree_t_m_avg{2,sd}-ttree_d_m_avg{2,sd}, smoothness)-smooth(dtree_t_m_avg{2,sd}-dtree_d_m_avg{2,sd}, smoothness));
                    diff_m = diff_m(t_relev);diff_m(diff_m<0)=0;
                    %                     i = find(diff_m<0, 1);diff_m(i:end)=0;
                    AUC_m_avg(sd) = sum(diff_m);
                    %             AUC_m(sd) = sum(smooth(t_connected_m_avg{sd}(t_relev)-t_unconnected_m_avg{sd}(t_relev), smoothness));
                end
                hold on; plot(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_m_avg, 'g')
            end
            %     legend({'no mask', 'masked'})
            xlabel('cue delay(ms)')
            title('Area Under the difference curve from stim end to stim end+300ms')
            e = std(AUC_m(:,good_chann),1,2)'/2;
            [hAx,hLine1,hLine2] = plotyy(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_m_avg,[Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)) 200],  [perf_m(2,1:min(numel(Poss_Cue_Time),6))./sum(perf_m(:,1:min(numel(Poss_Cue_Time),6))) perf(2)/sum(perf)]);
            hold on; errorbar(Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)), AUC_m_avg,e, 'g')
            hLine2.Color = [1 0 0];
            hLine1.Color = [0 .8 0];
            hAx(1).YLim = [0 80];
            hAx(2).YLim = [.5 .9];
            legend([h1;hLine1;hLine2],'no cue','cue','behavior','Location','northeast');
            
            filename = ['Diff_InfoTargNonTarg_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
            set(gcf, 'PaperPositionMode', 'auto')
            print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    
    
    neur_decay = [AUC_m_avg AUC_avg];
    cue_time = [Poss_Cue_Time(1:min(numel(Poss_Cue_Time),6)) NaN];
    behav_dacay = [perf_m(2,1:min(numel(Poss_Cue_Time),6))./sum(perf_m(:,1:min(numel(Poss_Cue_Time),6))) perf(2)/sum(perf)];
    filename = ['Darwin_varCue_NeurBehav_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    save(filename, 'neur_decay', 'behav_dacay', 'cue_time')
    %%
    %     compute integral of diff between connected and unconn curves over the
    %     interval 50-400
    t_relev = find(tb>.01 & tb<0.5);
    
    
    
    a =100:10:500;
    cum_int = zeros(1, numel(a));
    for t = 1:numel(a)
        t_relev = find(tb>.01 & tb<a(t)/1000);
        cum_int(t) = sum(ttree_t_avg{2}(t_relev)-dtree_t_avg{2}(t_relev));
    end
    figure; plot(a, cum_int)
    % 100ms of information :AUC
    t_relev = find(tb>0.1 & tb<.2);
    %     value = .5*(sum(t_connected_avg(t_relev)-t_unconnected_avg(t_relev))+ sum(t_connected_avg{6}(t_relev)-t_unconnected_avg{6}(t_relev))) ;%- sum(t_connected_avg{3}(t_relev)-t_unconnected_avg{3}(t_relev));
    
    %% all diff curves on same graph
    figure;hold on
    smoothness = 40;
    col = colormap('hot');
    col = col(1:10:end,:);
    for sd = 1:min(numel(Poss_Cue_Time),6)
        plot(tb, smooth(ttree_t_m_avg{2,sd}-dtree_t_m_avg{2,sd}, smoothness), 'Color', col(sd,:))
    end
    plot(tb, smooth(ttree_t_avg{2}-dtree_t_avg{2}, smoothness), 'k', 'LineWidth',2); xlim([-.3 1.2]);
    
    
    
else
    figure;
    for h = 1:2
        background_avg{h} = 0*background{25,h} ;n_background_avg(h)=0;
        ttree_t_avg{h} = 0*background{25,h};n_ttree_t_avg(h)=0;
        ttree_d_avg{h} = 0*background{25,h};n_ttree_d_avg(h)=0;
        dtree_t_avg{h} = 0*background{25,h};n_dtree_t_avg(h)=0;
        dtree_d_avg{h} = 0*background{25,h};n_dtree_d_avg(h)=0;
    end
    if plot_mask
        for h = 1:2
            background_m_avg{h} = 0*background{25,h} ;n_background_m_avg(h)=0;
            ttree_t_m_avg{h} = 0*background{25,h};n_ttree_t_m_avg(h)=0;
            ttree_d_m_avg{h} = 0*background{25,h};n_ttree_d_m_avg(h)=0;
            dtree_t_m_avg{h} = 0*background{25,h};n_dtree_t_m_avg(h)=0;
            dtree_d_m_avg{h} = 0*background{25,h};n_dtree_d_m_avg(h)=0;
        end
    end
    
    for n = good_chann
        
        for h = 1:2
            try background_avg{h} = nansum([background_avg{h} ,background{n,h}],2);end; if ~isnan(sum(background{n,h}));n_background_avg(h) = n_background_avg(h)+1;end
            
            try ttree_t_avg{h} = nansum([ttree_t_avg{h},ttree_t{n,h}],2);end; if ~isnan(sum(ttree_t{n,h}));n_ttree_t_avg=n_ttree_t_avg+1;end
            try ttree_d_avg{h} = nansum([ttree_d_avg{h} , ttree_d{n,h}],2);end; if ~isnan(sum(ttree_d{n,h}));n_ttree_d_avg=n_ttree_d_avg+1;end
            try dtree_t_avg{h} = nansum([dtree_t_avg{h}, dtree_t{n,h}],2);end; if ~isnan(sum(dtree_t{n,h}));n_dtree_t_avg=n_dtree_t_avg+1;end
            try dtree_d_avg{h} = nansum([dtree_d_avg{h} , dtree_d{n,h}],2);end; if ~isnan(sum(dtree_d{n,h}));n_dtree_d_avg=n_dtree_d_avg+1;end
        end
        if plot_mask
            for h = 1:2
                try background_m_avg{h} = nansum([background_m_avg{h} ,background_m{n,h}],2);end; if ~isnan(sum(background_m{n,h}));n_background_m_avg=n_background_m_avg+1;end
                try ttree_t_m_avg{h} = nansum([ttree_t_m_avg{h} ,ttree_t_m{n,h} ],2);end; if ~isnan(sum(ttree_t_m{n,h}));n_ttree_t_m_avg=n_ttree_t_m_avg+1;end
                try ttree_d_m_avg{h} = nansum([ttree_d_m_avg{h} , ttree_d_m{n,h}],2);end; if ~isnan(sum(ttree_d_m{n,h}));n_ttree_d_m_avg=n_ttree_d_m_avg+1;end
                try dtree_t_m_avg{h} = nansum([dtree_t_m_avg{h}, dtree_t_m{n,h}],2);end; if ~isnan(sum(dtree_t_m{n,h}));n_dtree_t_m_avg=n_dtree_t_m_avg+1;end
                try dtree_d_m_avg{h} = nansum([dtree_d_m_avg{h} , dtree_d_m{n,h}],2);end; if ~isnan(sum(dtree_d_m{n,h}));n_dtree_d_m_avg=n_dtree_d_m_avg+1;end
            end
        end
        
    end
    for h = 1:2
        background_avg{h} = background_avg{h}/n_background_avg(h);
        ttree_t_avg{h} = ttree_t_avg{h}/n_ttree_t_avg(h);
        ttree_d_avg{h} = ttree_d_avg{h}/n_ttree_d_avg(h);
        dtree_t_avg{h} = dtree_t_avg{h}/n_dtree_t_avg(h);
        dtree_d_avg{h} = dtree_d_avg{h}/n_dtree_d_avg(h);
    end
    if plot_mask
        
        for h = 1:2
            background_m_avg{h} = background_m_avg{h}/n_background_m_avg(h);
            ttree_t_m_avg{h} = ttree_t_m_avg{h}/n_ttree_t_m_avg(h);
            ttree_d_m_avg{h} = ttree_d_m_avg{h}/n_ttree_d_m_avg(h);
            dtree_t_m_avg{h} = dtree_t_m_avg{h}/n_dtree_t_m_avg(h);
            dtree_d_m_avg{h} = dtree_d_m_avg{h}/n_dtree_d_m_avg(h);
        end
    end
    
    if plot_mask; subplot(1,2,1);end; hold on
%     plot(tb,smooth(background_avg{1},smoothness),'k');
%     plot(tb,smooth(ttree_t_avg{1},smoothness),'b');
%     plot(tb,smooth(ttree_d_avg{1},smoothness),'r');
%     plot(tb,smooth(dtree_t_avg{1},smoothness),'c');
%     plot(tb,smooth(dtree_d_avg{1},smoothness),'m');
    legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'})
    plot(tb,smooth(background_avg{2},smoothness),'k', 'LineWidth',2);
    plot(tb,smooth(ttree_t_avg{2},smoothness),'b', 'LineWidth',2);
    plot(tb,smooth(ttree_d_avg{2},smoothness),'r', 'LineWidth',2);
    plot(tb,smooth(dtree_t_avg{2},smoothness),'c', 'LineWidth',2);
    plot(tb,smooth(dtree_d_avg{2},smoothness),'m', 'LineWidth',2); xlim([-.3 1.2]);
    if plot_mask
        subplot(1,2,2); hold on
%         plot(tb,smooth(background_m_avg{1},smoothness),'k');
%         plot(tb,smooth(ttree_t_m_avg{1},smoothness),'b');
%         plot(tb,smooth(ttree_d_m_avg{1},smoothness),'r');
%         plot(tb,smooth(dtree_t_m_avg{1},smoothness),'c');
%         plot(tb,smooth(dtree_d_m_avg{1},smoothness),'m');
        legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'})
        plot(tb,smooth(background_m_avg{2},smoothness),'k', 'LineWidth',2);
        plot(tb,smooth(ttree_t_m_avg{2},smoothness),'b', 'LineWidth',2);
        plot(tb,smooth(ttree_d_m_avg{2},smoothness),'r', 'LineWidth',2);
        plot(tb,smooth(dtree_t_m_avg{2},smoothness),'c', 'LineWidth',2);
        plot(tb,smooth(dtree_d_m_avg{2},smoothness),'m', 'LineWidth',2); xlim([-.3 1.2]);
        
    end
end


% figure,bar(SNR)
% figure,bar(peak)
%
% figure; for t = 1:155; plot(tb, e(:,t)); pause(0.1);end
