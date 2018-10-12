% function to analyse the redo of the decay of iconic memory experiment in
% Darwin - sessions 20180808 and later. based on
% figana_chanbychan_MontyIcoMemo_Decay_new.m
% function created by Catherine Wacongne (catherine.waco@gmail.com)
% final: sessions 136:150 were recorded with incorrect timings
% session 151 onwards has good timings

close all
clearvars
dbstop if error

plot_corr = 0; % if 1 analyses only correct trials. Otherwise takes all trials where the green fp was shown
integr_method = 'relu'; %option for computation of the integral ('abs' or 'relu'(<0 = 0), any other string gives regular integration) 
integ_dur = 0.50;
smoothness = 20;
baseline = [-.5 -.3];
x_limits = [-.5 1.2];
SNR_threshold = 1;

[extractdir, sessions, rawdir, SaveDir] = infoDir_DarwinIcoMemo;
info = Log_DarwinIcoMemo;


load([extractdir, sessions{end} '\EVT_', info(numel(sessions)).Tankname,'_block-',num2str(info(numel(sessions)).goodblocs(1))])
i = find(strcmp('ENV1',{EVENT.strms.name}));
SF = EVENT.strms(i).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);
StimCondStr = {'Short', 'Long'};
Cond = cell(96,2,2,2);
ci = zeros(96,2);

%132 sessions with fix stim dur and var mask onset
% timing issues fixed from session 151 onwards
Sessions = 151:160; %136:150; % 131:135; %132:135;%
SNR = zeros(1,48);

% target connected in RF: 6
% target unconnected in RF: 5 (target pair distractor)
curveInRF = 6; % pair in RF, curve in RF
curveOutRF = 5; % pair in RF, curve out RF

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
                        trials = trials(1:size(e,2));     
                    end
                    
                    %%% REMOVE THE FOLLOWING LINE IF THERE WAS NO DIODE
                    if isfield(LOG,'diode')
                        trials = trials(find(LOG.diode.goodtrials));
                        e = e(:,find(LOG.diode.goodtrials));
                    end
                    
                    days = [days, s*ones(1,size(e,2))];
                    if isempty(AllTrials)
                        AllTrials=e;
                        %                         connected = LOG.conn_codes(trials,1);
                        targ = LOG.targ_num(trials);
                        correct = LOG.correct(trials)+1;
                        cue_time = LOG.MaskDelay(trials);
                        stim_dur = LOG.Stim_dur(trials); % Stimulus duration
                        islong = LOG.conn_codes(trials,curveInRF)'-1;
                        angle = LOG.angle(trials);
                        if s>124
                            mask = LOG.isMask(trials);
                        else
                            error('script not valid for this session')
                        end
                    else
                        AllTrials = [AllTrials e];
                        %                         connected = [connected LOG.conn_codes(trials,1)];
                        targ = [targ LOG.targ_num(trials)];
                        correct = [correct LOG.correct(trials)+1];
                        cue_time = [cue_time LOG.MaskDelay(trials)];
                        stim_dur = [stim_dur LOG.Stim_dur(trials)]; % Stimulus duration
                        islong = [islong LOG.conn_codes(trials,curveInRF)'-1];
                        angle = [angle LOG.angle(trials)];
                        if s>124
                            mask = [mask LOG.isMask(trials)];
                        else
                            error('script not valid for this session')
                        end
                    end
                    
                    for  h = 1:2
                        pct = unique(cue_time);
                        for ct = 1:numel(pct)
                            ind = find(LOG.isMask(trials)==1 & (LOG.correct(trials)+1)==h  & LOG.MaskDelay(trials)==pct(ct));
                            perf(s, h, ct) = numel(ind);
                        end
                    end
                end
            end
            
            Poss_dur = unique(stim_dur);
            Poss_CueTime = unique(cue_time);
            Poss_Angle = unique(angle);
            bf = find(tb<baseline(2));
            if 0 %s<115
                pf = find(tb > 0.05-LOG.Precue_dur(1)/1e3 & tb < 0.1-LOG.Precue_dur(1)/1e3);
                pf2 = find(tb > 0.05-LOG.Precue_dur(1)/1e3 & tb < 0.09-LOG.Precue_dur(1)/1e3);
            else
                pf = find(tb > 0.05 & tb < 0.1);
                pf2 = find(tb > 0.05 & tb < 0.09);
            end
            
            hitNorm =(islong==1 & correct==2);
            [AllTrials2,eraseday, isgood] = NormalizeActiv(AllTrials, n, days,hitNorm ,bf, pf);
            aregood(:,n) = isgood;
            SNR(n) = nanmax(nanmean(AllTrials2(pf2,hitNorm==1),2))/nanmean(nanstd(AllTrials2(bf,:),0,2));
            AllTrials = AllTrials2;
           
            ind_a = cell(1,numel(Poss_Angle));
            for h = 1:2
                for ct = 1:numel(Poss_CueTime)
                    
%                     % average per condition: connected not masked trials
%                     for ang = 1:numel(Poss_Angle)
%                         if plot_corr
%                             ind_a{ang} = find(islong==1  & mask==0 & correct==h & cue_time==Poss_CueTime(ct) & angle==Poss_Angle(ang));%
%                         else
%                             ind_a{ang} = find(islong==1  & mask==0 & cue_time==Poss_CueTime(ct) & angle==Poss_Angle(ang));
%                         end
%                     end
%                     mintrialspercond = min(cellfun(@numel,ind_a));
%                     ind = [];
%                     for ang = 1:numel(Poss_Angle)
%                         rd = randperm(numel(ind_a{ang}));
%                         ind = [ind ind_a{ang}(rd(1:mintrialspercond))];
%                     end
%                     connected{n,h,ct} = nanmean(AllTrials(:,ind),2);
                    %                   
                    
                    
%                     % unconnected not masked trials
%                     for ang = 1:numel(Poss_Angle)
%                         if plot_corr
%                             ind_a{ang} = find(islong==0 & mask==0 & correct==h  & cue_time==Poss_CueTime(ct)& angle==Poss_Angle(ang));
%                         else
%                             ind_a{ang} = find(islong==0 & mask==0 & cue_time==Poss_CueTime(ct)& angle==Poss_Angle(ang));
%                         end
%                     end
%                     mintrialspercond = min(cellfun(@numel,ind_a));
%                     ind = [];
%                     for ang = 1:numel(Poss_Angle)
%                         rd = randperm(numel(ind_a{ang}));
%                         ind = [ind ind_a{ang}(rd(1:mintrialspercond))];
%                     end
%                     unconnected{n,h,ct} = nanmean(AllTrials(:,ind),2);
                    %                     
                    
                    
                    % connected  masked trials
                    for ang = 1:numel(Poss_Angle)
                        if plot_corr
                            ind_a{ang} = find(islong==1 & mask==1 & correct==h & cue_time==Poss_CueTime(ct)& angle==Poss_Angle(ang));
                        else
                            ind_a{ang} = find(islong==1 & mask==1 & cue_time==Poss_CueTime(ct)& angle==Poss_Angle(ang));
                        end
                    end
                    mintrialspercond = min(cellfun(@numel,ind_a));
                    ind = [];
                    for ang = 1:numel(Poss_Angle)
                        rd = randperm(numel(ind_a{ang}));
                        ind = [ind ind_a{ang}(rd(1:mintrialspercond))];
                    end
                    connected_m{n,h,ct} = nanmean(AllTrials(:,ind),2);
                    %                     
                    
                    % unconnected masked trials
                    for ang = 1:numel(Poss_Angle)
                        if plot_corr
                            ind_a{ang} = find(islong==0 & mask==1 & correct==h & cue_time==Poss_CueTime(ct)& angle==Poss_Angle(ang));
                        else
                            ind_a{ang} = find(islong==0 & mask==1 & cue_time==Poss_CueTime(ct)& angle==Poss_Angle(ang));
                        end
                    end
                    mintrialspercond = min(cellfun(@numel,ind_a));
                    ind = [];
                    for ang = 1:numel(Poss_Angle)
                        rd = randperm(numel(ind_a{ang}));
                        ind = [ind ind_a{ang}(rd(1:mintrialspercond))];
                    end
                    unconnected_m{n,h,ct} = nanmean(AllTrials(:,ind),2);
                    %                     
%                     
%                     ind = find(mask==0 & correct==h  & cue_time==Poss_CueTime(ct));
%                     perf(h, ct) = numel(ind);
                    ind = find(mask==1 & correct==h  & cue_time==Poss_CueTime(ct));
                    perf_m(h, ct) = numel(ind);
                    
                end
                
                
            end
            
        end
        
        %% plot the conditions for each channel
%         AUC_chan = zeros(48,numel(Poss_CueTime));
        AUC_mask_chan = zeros(48,numel(Poss_CueTime));
        close all
        
        
        f2 =figure('Position', [100 50 1400 900]);
        %n,h,sd,m
        
%         AUC = zeros(min(numel(Poss_CueTime),6),48);
        AUC_m = zeros(min(numel(Poss_CueTime),6),48);
        for n = ana_chan
            
            figure(f2);clf
            for ct = 1:min(numel(Poss_CueTime),6)
                
                figure(f2)
                subplot(2,4,ct);hold on; title([num2str(Poss_CueTime(ct)), ' ms, chan ' num2str(n)])
%                 plot(tb, smooth(connected{n,2,ct}, smoothness), 'b', 'LineWidth',1); 
%                 plot(tb, smooth(unconnected{n,2,ct}, smoothness), 'r', 'LineWidth',1);
                plot(tb, smooth(connected_m{n,2,ct}, smoothness), 'b', 'LineWidth',1);
                try  ylim([min(smooth(connected_m{n,2,ct},smoothness))-(max(smooth(connected_m{n,2,ct},smoothness))-min(smooth(connected_m{n,2,ct},smoothness))),  max(smooth(connected_m{n,2,ct},smoothness))+.2*(max(smooth(connected_m{n,2,ct},smoothness))-min(smooth(connected_m{n,2,ct},smoothness)))]);end
                [hAx,hLine1,hLine2] = plotyy(tb, smooth(unconnected_m{n,2,ct}, smoothness), tb, [smooth(connected_m{n,2,ct}-unconnected_m{n,2,ct}, smoothness)']);
                diffcu =smooth(connected_m{n,2,ct}-unconnected_m{n,2,ct},smoothness);
                hAx(2).XLim = x_limits;
                hAx(1).XLim = x_limits;
                hLine1.LineWidth = 1;
                hLine1.Color = [1 0 0];
                try hAx(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
                hLine2(1).Color=[0 1 0];   
                hLine2(1).LineWidth = 1;
                
                
                
                axes(hAx(2));line(hAx(2).XLim,[0 0]);
                p=patch([0.1 0.1+integ_dur 0.1+integ_dur 0.1],[-10 -10 20 20],'r');
                set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
                
            end
            
            
            
            AUC_m       = zeros(1,numel(Poss_CueTime));
            OffResp     = zeros(1,numel(Poss_CueTime));
            UncoOnResp  = zeros(1,numel(Poss_CueTime));
            t_relev = find(tb>(0.01) & tb<(0.01+integ_dur));
            for ct = 1:numel(Poss_CueTime)

                
                diff =smooth(connected_m{n,2,ct}, smoothness)-smooth(unconnected_m{n,2,ct}, smoothness);
                switch integr_method
                    case 'abs'
                        diff = abs(diff);
                    case 'relu'
                       diff(diff<0) = 0;
                end
                diff = diff(t_relev);%diff(diff<0)=0;
                AUC_m(ct) = sum(diff);
                
%                 tbl = find(tb>((1e-3*Poss_CueTime(ct))+0.01) & tb<((1e-3*Poss_CueTime(ct)) +0.03));
%                 blact = smooth(connected_m{n,2,ct}(tbl)-unconnected_m{n,2,ct}(tbl),20);
%                 pre_off = mean(blact);
%                 pfskel = find(tb > 0.05-LOG.Precue_dur(1)/1e3 & tb < 0.12-LOG.Precue_dur(1)/1e3);
%                 UncoOnResp(ct) = max(smooth(unconnected_m{n,2,ct}(pfskel), 20));%
%                 tpeak = find(tb>((1e-3*Poss_CueTime(ct))+0.03) & tb<((1e-3*Poss_CueTime(ct)) +0.09));
%                 pact = smooth(connected_m{n,2,ct}(tpeak)-unconnected_m{n,2,ct}(tpeak),20);
%                 peakvalue = max(pact);
%                 OffResp(ct) = peakvalue - pre_off;
%                 OffRatio(ct) = peakvalue/pre_off;
                
            end
            subplot(2,4,8);hold on; plot(Poss_CueTime, AUC_m, 'Color', [0 .9 0]);
            
%             try
%                 flin = @(param,xval) param(1)+param(2)*xval;
%                 flin_m = @(param,xval) param(3)+param(2)*xval;
%                 [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_CueTime, Poss_CueTime}, {AUC,AUC_m}, {flin, flin_m}, [20 20 0.2]);
%                 x2 = 0:10:Poss_CueTime(end);
%                 y2 = beta(1)+beta(2)*x2;
%                 y2_m =  beta(3)+beta(2)*x2;
%                 plot(x2, y2, 'Color','k')
%                 plot(x2, y2_m, 'Color','g')
%                 [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked', 'Location', 'NorthWest');
%                 LegendPos = get(legend_h,'Position');
%                 EstimWorthNeuronal_chann(n) = (beta(1)-beta(3))/beta(2);
%                 dim = [LegendPos(1) LegendPos(2)-LegendPos(4) .3 LegendPos(4)/2];
%                 str = ['Estimated worth: ' num2str(EstimWorthNeuronal_chann(n),3) 'ms'];
%                 annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
%                 
%             catch
%                 EstimWorthNeuronal_chann(n) = NaN;
%             end
            
            
            
            AUC_mask_chan(n,:) = AUC_m;
            xlabel('stimulus duration')
            pause(0.02)
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-dpng',[SaveDir 'ConnUnconn_Chan' num2str(n) '_s' num2str(days(1)) '_' num2str(days(end)) '.png'])  
            
        end
        
        %% accross chan analysis 
        good_chann = find(SNR>SNR_threshold);
       
%             figure;plot(OffResp_chan(good_chann), EstimWorthNeuronal_chann(good_chann), 'o')
%             [R,P] = corrcoef(OffResp_chan(good_chann),EstimWorthNeuronal_chann(good_chann));
%             linear_fit = polyfit(OffResp_chan(good_chann),EstimWorthNeuronal_chann(good_chann),1);
%             x_fit = 0:0.1:(1.1*max(OffResp_chan(good_chann)));
%             hold on; plot(x_fit, linear_fit(2)+linear_fit(1)*x_fit);
%             xlabel('Amplitude of Offset Response')
%             ylabel('Estimated Worth')
%             title({'correlation between offset response and', 'neuronal worth at single channel level'})
%             filename = ['CorrOffsetWorth_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
%             set(gcf, 'PaperPositionMode', 'auto')
%             print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
%          
%         
%         
%             figure;plot(UncoOnResp_chan(good_chann), EstimWorthNeuronal_chann(good_chann), 'o')
%             [R,P] = corrcoef(UncoOnResp_chan(good_chann),EstimWorthNeuronal_chann(good_chann));
%             linear_fit = polyfit(UncoOnResp_chan(good_chann),EstimWorthNeuronal_chann(good_chann),1);
%             x_fit = 0:0.1:(1.1*max(UncoOnResp_chan(good_chann)));
%             hold on; plot(x_fit, linear_fit(2)+linear_fit(1)*x_fit);
%             xlabel('Amplitude of Unconnected Onset Resp')
%             ylabel('Estimated Worth')
%             title({'correlation between unconnected onset response and', 'neuronal worth at single channel level'})
%             filename = ['CorrOffsetWorth_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
%             set(gcf, 'PaperPositionMode', 'auto')
%             print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
%      
%         
%         
%        
%             figure;plot(UncoOnResp_chan(good_chann), OffResp_chan(good_chann), 'o')
%             [R,P] = corrcoef(UncoOnResp_chan(good_chann),OffResp_chan(good_chann));
%             linear_fit = polyfit(UncoOnResp_chan(good_chann),OffResp_chan(good_chann),1);
%             x_fit = 0:0.1:(1.1*max(UncoOnResp_chan(good_chann)));
%             hold on; plot(x_fit, linear_fit(2)+linear_fit(1)*x_fit);
%             xlabel('Amplitude of Unconnected Onset Resp')
%             ylabel('Amplitude of Offset Response')
%             title({'correlation between unconnected onset response and', 'offset response at single channel level'})
%             filename = ['CorrOffsetWorth_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
%             set(gcf, 'PaperPositionMode', 'auto')
%             print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
%         
        
        %% plot performances
        good_chann = find(SNR>SNR_threshold);
        
        figure; hold on
        title('Behavioral performance on the iconic memory task');

%         hold on;  plot(Poss_CueTime(1:min(numel(Poss_CueTime),6)),  perf(2,1:min(numel(Poss_CueTime),6))./sum(perf(:,1:min(numel(Poss_CueTime),6))),'b')
        plot(Poss_CueTime(1:min(numel(Poss_CueTime),6)),  perf_m(2,1:min(numel(Poss_CueTime),6))./sum(perf_m(:,1:min(numel(Poss_CueTime),6))),'r')
        hold all
        for b=1:size(perf,1)
            plot(Poss_CueTime(1:min(numel(Poss_CueTime),6)), squeeze(perf(b,2,1:min(numel(Poss_CueTime),6))./sum(perf(b,:,1:min(numel(Poss_CueTime),6)))),'k')
        end
%         fsigm = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(2)-xval)*param(3)));
%         fsigm_m = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(4)-xval)*param(3)));
%         [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_CueTime, Poss_CueTime}, {perf(2,:)./sum(perf),perf_m(2,:)./sum(perf_m)}, {fsigm, fsigm_m}, [.8  150 0.005 0.005]);

%         [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked', 'Location', 'SouthEast');
%         LegendPos = get(legend_h,'Position');
%         set(legend_h, 'Location', 'NorthWest')
%         EstimWorthBehav = (beta(4)-beta(2));
%         dim = [LegendPos(1)-LegendPos(3) LegendPos(2) .3 LegendPos(4)/2];
%         str = ['Estimated worth: ' num2str(EstimWorthBehav,3) 'ms'];
%         annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
        
        
%         x1 = 0:10:Poss_CueTime(end);
%         y1 = 0.5+(beta(1)-0.5)./(1+10.^((beta(2)-x1)*beta(3)));
%         y1_m = 0.5+(beta(1)-0.5)./(1+10.^((beta(4)-x1)*beta(3)));
%         plot(x1, y1, 'Color','b')
%         plot(x1, y1_m, 'Color','r')
        xlabel('mask delay (ms)')

        filename = ['performances_' StimCondStr{stimcond} 'Stim_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        
        %% plot the average over an array for all the conditions
        f2 = figure('Position', [100 50 1400 900]);f3 = figure('Position', [100 50 1400 500]);
        for ct = 1:min(numel(Poss_CueTime),6)
            % compute the average
            for h = 1:2
%                 connected_avg{h,ct} = 0*connected{ana_chan(1),1};n_connected_avg(h,ct)=0;
%                 unconnected_avg{h,ct} = 0*connected{ana_chan(1),1};n_unconnected_avg(h,ct)=0;
                connected_m_avg{h,ct} = 0*connected_m{ana_chan(1),1};n_connected_m_avg(h,ct)=0;
                unconnected_m_avg{h,ct} = 0*connected_m{ana_chan(1),1};n_unconnected_m_avg(h,ct)=0;
            end
            for n = good_chann
                for h = 1:2
%                     try  connected_avg{h,ct} =  nansum([connected_avg{h,ct} ,connected{n,h,ct}],2);end;if ~isnan(sum(connected{n,h,ct}));n_connected_avg(h,ct) = n_connected_avg(h,ct)+1;end
%                     try unconnected_avg{h,ct} =  nansum([unconnected_avg{h,ct} ,unconnected{n,h,ct} ],2);end;if ~isnan(sum(unconnected{n,h,ct}));n_unconnected_avg(h,ct) = n_unconnected_avg(h,ct)+1;end
                    
                    try connected_m_avg{h,ct} =  nansum([connected_m_avg{h,ct},connected_m{n,h,ct}],2);end;if ~isnan(sum(connected_m{n,h,ct}));n_connected_m_avg(h,ct) = n_connected_m_avg(h,ct)+1;end
                    try unconnected_m_avg{h,ct} =  nansum([unconnected_m_avg{h,ct},unconnected_m{n,h,ct}],2);end;if ~isnan(sum(unconnected_m{n,h,ct}));n_unconnected_m_avg(h,ct) = n_unconnected_m_avg(h,ct)+1;end
                end
                
            end
            for h = 1:2
%                 connected_avg{h,ct} = connected_avg{h,ct}/n_connected_avg(h,ct);
%                 unconnected_avg{h,ct} = unconnected_avg{h,ct}/n_unconnected_avg(h,ct);
                connected_m_avg{h,ct} = connected_m_avg{h,ct}/n_connected_m_avg(h,ct);
                unconnected_m_avg{h,ct} = unconnected_m_avg{h,ct}/n_unconnected_m_avg(h,ct);
            end
            
            % plot
            figure(f2);subplot(2,3,ct);hold on; title(['cued at ' num2str(Poss_CueTime(ct)), ' ms'])%
%             plot(tb,smooth(connected_avg{1,ct},smoothness),'b');
%             plot(tb,smooth(unconnected_avg{1,ct},smoothness),'r');
            plot(tb,smooth(connected_m_avg{1,ct},smoothness),'b');
            plot(tb,smooth(unconnected_m_avg{1,ct},smoothness),'r');
            if ct == 1; legend({'target', 'target tree unconnected', 'distractor', 'distractor unconnected'});end
%             plot(tb,smooth(connected_avg{2,ct},smoothness),'b', 'LineWidth',2);
%             plot(tb,smooth(unconnected_avg{2,ct},smoothness),'r', 'LineWidth',2);
            plot(tb,smooth(connected_m_avg{2,ct},smoothness),'b', 'LineWidth',2);
            plot(tb,smooth(unconnected_m_avg{2,ct},smoothness),'r', 'LineWidth',2); xlim([-.4 0.9]);
            ylim([-.2 1.2])
        end
        
        figure(f2)
        filename = ['NeurActiv_withStimDur_Cue_' StimCondStr{stimcond} 'Stim_s' num2str(Sessions(1)) '_' num2str(Sessions(end)) '_chan' num2str(ana_chan(1)) '_' num2str(ana_chan(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        
        %% plot all the cue times for each conditions
        figure(f3);clf
%         subplot(2,2,1);hold on; title(['Connected'])%
%         labels = cell(1,numel(Poss_CueTime));
%         for ct = 1:numel(Poss_CueTime)
%             plot(tb,smooth(connected_avg{2,ct},smoothness),'Color', [max(0,.2*ct-.4) max(0,.2*ct-.4) min(1,.4*ct)]);
%             labels{ct} = ['Cue at ' num2str(Poss_CueTime(ct)) 'ms'];
%         end
%         xlim([-.2 .9]);
%         legend(labels)
%         subplot(2,2,2);hold on; title(['Unconnected'])%
%         for ct = 1:numel(Poss_CueTime)
%             plot(tb,smooth(unconnected_avg{2,ct},smoothness),'Color', [min(1,.4*ct) max(0,.2*ct-.4) max(0,.2*ct-.4)]);
%         end
%         xlim([-.2 .9]);
        subplot(1,2,1);hold on; title(['Connected masked'])%
        for ct = 1:numel(Poss_CueTime)
            plot(tb,smooth(connected_m_avg{2,ct},smoothness),'Color', [max(0,.2*ct-.4) max(0,.2*ct-.4) min(1,.4*ct)]);
        end
        xlim([-.2 .9]);
        subplot(1,2,2);hold on; title(['Unconnected masked'])%
        for ct = 1:numel(Poss_CueTime)
            plot(tb,smooth(unconnected_m_avg{2,ct},smoothness),'Color', [min(1,.4*ct) max(0,.2*ct-.4) max(0,.2*ct-.4)]);
        end
        xlim([-.2 .9]);
        filename = ['NeurActiv_withCond_Cue_' 's' num2str(Sessions(1)) '_' num2str(Sessions(end)) '_chan' num2str(ana_chan(1)) '_' num2str(ana_chan(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        
       
        %% plot the difference between connected and unconnected for target and distr pair
        figure('Position', [100 50 1400 900])
%         AUC_avg = zeros(1,min(numel(Poss_CueTime),6));
        AUC_m_avg = zeros(1,min(numel(Poss_CueTime),6));
        
        for cd = 1:min(numel(Poss_CueTime),7)
            topplot = cd+4*floor((cd-1)/4);
            bottomplot = cd+4*floor((cd-1)/4)+4;
           
            subplot(4,4,topplot);hold on; title(['Mask onset: ' num2str(Poss_CueTime(cd)), ' ms']); xlabel('time(s)'); ylabel('normalized activity')
            ylim1 = [-2 2];
            if 1
                a = plot(tb, smooth(connected_m_avg{2,cd}, smoothness), 'b', 'LineWidth',1); ylim(ylim1); hold all;
                b = plot(tb, smooth(unconnected_m_avg{2,cd}, smoothness),'r');
                
                subplot(4,4,bottomplot);
                plot(tb, smooth(connected_m_avg{2,cd}-unconnected_m_avg{2,cd}, smoothness),'Color',[0 1 0]); hold all;
                xlim([-0.2 0.9]);
                ylim([-0.2 1]);
                
                
            end
            
            %set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
%             if cd==1
%                 legend([b, hLine12, hLine22], {'no mask connected', 'no mask unconnected', 'conn - unconn'})
%             end
            diff = smooth(connected_m_avg{2,cd}-unconnected_m_avg{2,cd}, smoothness);
            switch integr_method
                case 'abs'
                    diff = abs(diff);
                case 'relu'
                    diff(diff<0) = 0;
            end
            difftotal = diff;
            diff = diff(t_relev);%diff(diff<0)=0;
            
            trelevnew = find(diff>0);
            
            qx = tb(t_relev);
            qy = diff;
            q1 = find(diff~=0,1,'Last');
            %area(qx(1:q1),qy(1:q1));
            
            
            greenarea = area(tb(t_relev),diff,0);
            
            %greenarea = area(tb(t_relev(trelevnew)),difftotal(trelevnew),0);
            set(greenarea,'FaceColor',[0,1,0],'FaceAlpha',0.5,'EdgeColor','none'); hold on;

            xlabel('time(s)');
            ylabel('normalized activity');
            
            AUC_m_avg(cd) = sum(diff);
            
            
            stimend = Poss_dur/1000;
            maskon = stimend+Poss_CueTime(cd)/1000;
            maskdur = 100/1000;
            xlim([-0.2 0.9]);        
            ylim([-0.2 1]);
            subplot(4,4,topplot)
            
            pstim=patch([0 0+stimend 0+stimend 0],[-10 -10 20 20],'b');
            set(pstim,'FaceAlpha',0.2, 'EdgeColor', 'none' );
            pstim=patch([maskon maskon+maskdur maskon+maskdur maskon],[-10 -10 20 20],'r');
            set(pstim,'FaceAlpha',0.2, 'EdgeColor', 'none' )
            xlim([-0.2 0.9]);
            ylim([-0.2 1]);

        end
        
        
        
        subplot(2,4,8);hold on;
%         hold on; plot(Poss_CueTime(1:min(numel(Poss_CueTime),6)), AUC_avg, 'ok'); plot(Poss_CueTime(1:min(numel(Poss_CueTime),6)), AUC_m_avg, 'og')
        
        xlabel('mask delay(ms)')
        title('Area Under the difference curve from stim end to fp green')
%         e = std(AUC_chan(good_chann,:),1,1)'/2;
        e_m = std(AUC_mask_chan(good_chann,:),1,1)'/2;
        
        perf_m_average = [perf_m(2,1:min(numel(Poss_CueTime),6))./sum(perf_m(:,1:min(numel(Poss_CueTime),6))) ];
        perf_m_range = max(perf_m_average)-min(perf_m_average); 
        AUC_m_range = max(AUC_m_avg)-min(AUC_m_avg);
        [hAx,hLine1,hLine2] = plotyy(Poss_CueTime(1:min(numel(Poss_CueTime),6)), AUC_m_avg,[Poss_CueTime(1:min(numel(Poss_CueTime),6)) ],  perf_m_average);
       
        set(hAx(1),'YLim',[min(AUC_m_avg)-0.5*AUC_m_range max(AUC_m_avg)+0.5*AUC_m_range]);
        set(hAx(2),'YLim',[min(perf_m_average)-0.5*perf_m_range max(perf_m_average)+0.5*perf_m_range]);
        
        hold on; errorbar(Poss_CueTime(1:min(numel(Poss_CueTime),6)), AUC_m_avg,e_m, '-og') %errorbar(Poss_CueTime(1:min(numel(Poss_CueTime),6)), AUC_avg,e, 'ok');
        %         hLine2.Color = [1 0 0];
        %         hLine1.Color = [0 .8 0];
        %         hAx(1).YLim = [0 160];
        %         hAx(2).YLim = [.5 .8];
        %         legend([hLine1;hLine2],'neuronal','behavior','Location','northeast');
        
        % --- STATISTICS
        figure
        subplot(1,2,1); 
        plot(Poss_CueTime,perf_m_average);
        subplot(1,2,2); 
        plot(Poss_CueTime,AUC_m_avg);
        
        % linear regression
        [b,bint,r,rint,stats] = regress(AUC_m_avg',[ones(length(perf_m_average),1) perf_m_average']);
        
        %% behavior
%         plot(Poss_CueTime(1:min(numel(Poss_CueTime),6)),  perf_m(2,1:min(numel(Poss_CueTime),6))./sum(perf_m(:,1:min(numel(Poss_CueTime),6))),'r')
%         [hAx,hLine1,hLine2] = plotyy(Poss_CueTime(1:min(numel(Poss_CueTime),6)), AUC_avg,[Poss_CueTime(1:min(numel(Poss_CueTime),6)) ],  [perf(2,1:min(numel(Poss_CueTime),6))./sum(perf(:,1:min(numel(Poss_CueTime),6))) ]);
        
%         flin = @(param,xval) param(1)+param(2)*xval;
%         flin_m = @(param,xval) param(3)+param(2)*xval;
%         [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_CueTime, Poss_CueTime}, {AUC_avg,AUC_m_avg}, {flin, flin_m}, [20 20 0.2]);
%         x2 = 0:10:Poss_CueTime(end);
%         y2 = beta(1)+beta(2)*x2;
%         y2_m =  beta(3)+beta(2)*x2;
%         plot(x2, y2, 'Color','k')
%         plot(x2, y2_m, 'Color','g')
%         [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked', 'Location', 'NorthWest');
%         LegendPos = get(legend_h,'Position');
%         EstimWorthNeuronal_chann(n) = (beta(1)-beta(3))/beta(2);
%         dim = [LegendPos(1) LegendPos(2)-LegendPos(4) .3 LegendPos(4)/2];
%         str = ['Estimated worth: ' num2str(EstimWorthNeuronal_chann(n),3) 'ms'];
%         annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
%                 
        
        
        
        
        filename = ['Diff_InfoTargNonTarg_' StimCondStr{stimcond} 'Stim_s' num2str(Sessions(1)) '_' num2str(Sessions(end)) '_chan' num2str(ana_chan(1)) '_' num2str(ana_chan(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
        
        
        neur_decay = [ AUC_m_avg];
        stim_durs = [Poss_CueTime(1:min(numel(Poss_CueTime),6)) NaN];
        behav_dacay = [perf_m(2,1:min(numel(Poss_CueTime),6))./sum(perf_m(:,1:min(numel(Poss_CueTime),6)))];
        filename = ['Monty_varCue_NeurBehav_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
        save(filename, 'neur_decay', 'behav_dacay', 'stim_durs')
        
        
       
    end
end


