% analysis function for the iconic memory experiment - control with 3 mask
% conditions at 0, 1frame and no mask and variable stimulus duration
% This script can be used for the analysis of the 
% (s29-30) 
% written by Catherine Wacongne 


close all
clear all
dbstop if error
plot_stim_dur = 1;
plot_mask = 1;
plot_chann= 1;
integ_dur = 0.3;
SNR_threshold = 1;
%directory of extracted data and stim log data
[extractdir, sessions, rawdir, SaveDir] = infoDir_MontyIcoMemo;
info = Log_MontyIcoMemo;
load([extractdir, sessions{1} '\EVT_', info(1).Tankname,'_block-',num2str(info(1).goodblocs(1))])


SF = EVENT.strms(1).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);


Cond = cell(48,2,2,2);
Sessions = 29:30;
for n = 1:48
    % Extract the trials from all blocs and conditions for each channel
    AllTrials = [];
    AllCond = [];
    days =[];
    for s = Sessions
        
        for bloc = 1:numel(info(s).goodblocs)
            %Read in the extracted data form the mapped channel
            load([extractdir sessions{s},'\Xtract_',info(s).Tankname,'_Block-',num2str(info(s).goodblocs(bloc)),'_',num2str(n)]);
            e = Env{1};
            load([rawdir,'ICOmemo_' info(s).Tankname '_B',num2str(info(s).goodblocs(bloc))])
            trials =find(LOG.target_presented>0);
            %             e = e(:,trials);
            equal = (size(e,2) == length(trials));
            if size(e,2)<length(trials)
                trials = trials(1:size(e,2));
            elseif size(e,2)>length(trials)
                disp(['skipping bloc ' num2str(bloc) ' in session ' num2str(s)])
                continue
            end
            
            isdrawn4_b = zeros(1,numel(trials));
            connected_b = zeros(1,numel(trials));
            for i =1: numel(trials)
                a = find(LOG.drawn_skeleton{trials(i)}==5);
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
                Hit = LOG.Hit(trials);
                stim_dur = LOG.Stim_dur(trials);
                Poss_dur = unique(stim_dur);
                mask = 2;
                if s>=29
                    mask = LOG.Masked(trials)+((LOG.Masked(trials).*LOG.Mask_delay(trials))>0)+1;
                else
                    mask = mask.*(LOG.Masked(trials)) +1;
                end
                
                %                 AllCond = MAT;
            else
                AllTrials = [AllTrials e];
                %                 AllCond = [AllCond;MAT];
                isdrawn4 = [isdrawn4 isdrawn4_b];
                connected = [connected connected_b];
                targ = [targ LOG.targ_num(trials)];
                Hit = [Hit LOG.Hit(trials)];
                stim_dur = [stim_dur LOG.Stim_dur(trials)];
                Poss_dur = [Poss_dur unique(stim_dur)];
                mask_d = 2;
                if s>=29
                    mask = [mask LOG.Masked(trials)+((LOG.Masked(trials).*LOG.Mask_delay(trials))>0)+1];
                else
                    mask =  [mask  mask_d.*(LOG.Masked(trials))+1];
                end
            end
            
        end
    end
    Poss_dur = unique(Poss_dur);
    % do some cleaning
    ev = reshape(AllTrials,size(AllTrials,1)*size(AllTrials,2),1);
    ev(abs(zscore(ev))>5) = NaN;
    AllTrials = reshape(ev,size(AllTrials,1),size(AllTrials,2));% figure;subplot(1,2,1);imagesc(AllTrials);
    
    %Z-score the trials
    T = nanmean(AllTrials);
    f1 = find(abs(zscore(T))>3);
    f2 = find(max(abs(AllTrials))>0.1);
    AllTrials(:,[f1,f2]) = NaN(size(AllTrials,1),length(f1)+length(f2));%subplot(1,2,2);imagesc(AllTrials);pause
    %     keyboard
    bf = find(tb<0);
    pf = find(tb > 0.03 & tb < 0.1);
    pf2 = find(tb > 0.05 & tb < 0.09);
    hitNorm =(targ==10 & Hit==2);
    %      AllTrials = NormalizeMUA(AllTrials, n, days,hitNorm ,bf, pf);
    AllTrials2 = NormalizeActiv(AllTrials, n, days,hitNorm ,bf, pf2);
    %      AllTrials = NormalizeMUA(AllTrials, n, days,hitNorm ,bf, pf);
    SNR(n) = nanmax(nanmean(AllTrials2(pf2,hitNorm==1),2))/nanmean(nanstd(AllTrials2(bf,:)));
    %     figure; subplot(1,2,1);imagesc(AllTrials);subplot(1,2,2); imagesc(AllTrials2);pause(0.1)
    AllTrials = AllTrials2;
    
    for h = 1:2
        if plot_stim_dur
            for sd = 1:numel(Poss_dur)
                %                             for m = 1:3
                % average per condition
                ind = find(isdrawn4==0 & Hit==h & stim_dur == Poss_dur(sd) & mask ==1 );%& stim_dur == Poss_dur(sd) & mask == m
                background{n,h,sd}  = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(targ==10 & Hit==h & stim_dur == Poss_dur(sd) & mask ==1);
                ttree_t{n,h,sd} = nanmean(AllTrials(:,ind),2);
                
                
                
                ind = find(targ==9 & Hit==h & stim_dur == Poss_dur(sd) & mask ==1);
                ttree_d{n,h,sd} = nanmean(AllTrials(:,ind),2);
                
                ind = find(connected==1 & targ~=9 & Hit==h & stim_dur == Poss_dur(sd)& mask ==1 );
                dtree_t{n,h,sd} = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(connected==2 & targ~=10 & Hit==h & stim_dur == Poss_dur(sd)& mask ==1 );
                dtree_d{n,h,sd} = nanmean(AllTrials(:,ind),2);
                
                ind = find( Hit==h & stim_dur == Poss_dur(sd) & mask ==1);
                perf(h,sd) = numel(ind);
                
                
                ind = find(connected==1 & stim_dur == Poss_dur(sd) & mask ==1);% & Hit==2);%
                t_connected{n,sd} =  nanmean(AllTrials(:,ind),2);
                ind = find(connected==2 & stim_dur == Poss_dur(sd) & mask ==1);%& Hit==2);%
                t_unconnected{n,sd} = nanmean(AllTrials(:,ind),2);
                if plot_mask
                    ind = find(isdrawn4==0 & Hit==h & stim_dur == Poss_dur(sd) & mask ==3 );%& stim_dur == Poss_dur(sd) & mask == m
                    background_m{n,h,sd}  = nanmean(AllTrials(:,ind),2);
                    
                    
                    ind = find(targ==10 & Hit==h & stim_dur == Poss_dur(sd) & mask ==3);
                    ttree_t_m{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    
                    
                    ind = find(targ==9 & Hit==h & stim_dur == Poss_dur(sd) & mask ==3);
                    ttree_d_m{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    ind = find(connected==1 & targ~=9 & Hit==h & stim_dur == Poss_dur(sd)& mask ==3 );
                    dtree_t_m{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    
                    ind = find(connected==2 & targ~=10 & Hit==h & stim_dur == Poss_dur(sd)& mask ==3 );
                    dtree_d_m{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    ind = find( Hit==h & stim_dur == Poss_dur(sd) & mask ==3);
                    perf_m(h,sd) = numel(ind);
                    
                    
                    ind = find(connected==1 & stim_dur == Poss_dur(sd) & mask ==3);% & Hit==2);%
                    t_connected_m{n,sd} =  nanmean(AllTrials(:,ind),2);
                    ind = find(connected==2 & stim_dur == Poss_dur(sd) & mask ==3);% & Hit==2);%
                    t_unconnected_m{n,sd} = nanmean(AllTrials(:,ind),2);
                    
                    
                    
                    % no delay mask
                    ind = find(isdrawn4==0 & Hit==h & stim_dur == Poss_dur(sd) & mask ==2 );%& stim_dur == Poss_dur(sd) & mask == m
                    background_m2{n,h,sd}  = nanmean(AllTrials(:,ind),2);
                    
                    
                    ind = find(targ==10 & Hit==h & stim_dur == Poss_dur(sd) & mask ==2);
                    ttree_t_m2{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    
                    
                    ind = find(targ==9 & Hit==h & stim_dur == Poss_dur(sd) & mask ==2);
                    ttree_d_m2{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    ind = find(connected==1 & targ~=9 & Hit==h & stim_dur == Poss_dur(sd)& mask ==2 );
                    dtree_t_m2{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    
                    ind = find(connected==2 & targ~=10 & Hit==h & stim_dur == Poss_dur(sd)& mask ==2 );
                    dtree_d_m2{n,h,sd} = nanmean(AllTrials(:,ind),2);
                    
                    ind = find( Hit==h & stim_dur == Poss_dur(sd) & mask ==2);
                    perf_m2(h,sd) = numel(ind);
                    
                    
                    ind = find(connected==1 & stim_dur == Poss_dur(sd) & mask ==2);% & Hit==2);%
                    t_connected_m2{n,sd} =  nanmean(AllTrials(:,ind),2);
                    ind = find(connected==2 & stim_dur == Poss_dur(sd) & mask ==2);% & Hit==2);%
                    t_unconnected_m2{n,sd} = nanmean(AllTrials(:,ind),2);
                    
                end
                
            end
        else
            ind = find(isdrawn4==0 & Hit==h & mask ==1 );%& stim_dur == Poss_dur(sd) & mask == m
            background{n,h}  = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(targ==10 & Hit==h & mask ==1 );
            ttree_t{n,h} = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(targ==9 & Hit==h & mask ==1 );
            ttree_d{n,h} = nanmean(AllTrials(:,ind),2);
            
            ind = find(connected==1 & targ~=9 & Hit==h & mask ==1);
            dtree_t{n,h} = nanmean(AllTrials(:,ind),2);
            
            
            ind = find(connected==2 & targ~=10 & Hit==h & mask ==1);
            dtree_d{n,h} = nanmean(AllTrials(:,ind),2);
            
            ind = find(connected==1   &  mask ==1);% & Hit==2);%
            t_connected{n} =  nanmean(AllTrials(:,ind),2);
            ind = find(connected==2   &  mask ==1);%& Hit==2);%
            t_unconnected{n} = nanmean(AllTrials(:,ind),2);
            
            if plot_mask
                ind = find(isdrawn4==0 & Hit==h & mask ==3 );%& stim_dur == Poss_dur(sd) & mask == m
                background_m{n,h}  = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(targ==10 & Hit==h & mask ==3 );
                ttree_t_m{n,h} = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(targ==9 & Hit==h & mask ==3 );
                ttree_d_m{n,h} = nanmean(AllTrials(:,ind),2);
                
                ind = find(connected==1 & targ~=9 & Hit==h & mask ==3);
                dtree_t_m{n,h} = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(connected==2 & targ~=10 & Hit==h & mask ==3);
                dtree_d_m{n,h} = nanmean(AllTrials(:,ind),2);
                
                ind = find(connected==1  & mask ==3);% & Hit==2);%
                t_connected_m{n} =  nanmean(AllTrials(:,ind),2);
                ind = find(connected==2  & mask ==3);% & Hit==2);%
                t_unconnected_m{n} = nanmean(AllTrials(:,ind),2);
                
            end
            
        end
    end
   
end

%% plot the conditions

close all
if plot_chann
    %     f1 =figure('Position', [100 50 1400 900]);
    f2 =figure('Position', [100 50 1400 900]);
    if ~plot_stim_dur
        f10 = figure('Position', [100 50 1400 900]);
    end
    %n,h,sd,m
    for n = 1:48
        
        
        if plot_stim_dur
            
            %             figure(f1);clf
            figure(f2);clf
            for sd = 1:min(numel(Poss_dur),7)
                
                ylim1 = [-2 1.5];
                
                figure(f2)
                subplot(2,4,sd);hold on; title([num2str(Poss_dur(sd)), ' ms, chan ' num2str(n)])
                
                if plot_mask
                    plot(tb, smooth(t_connected_m{n,sd}, 20), 'b', 'LineWidth',1); ylim(ylim1)%try  ylim([min(smooth(t_connected_m{n,sd},smoothness))-(max(smooth(t_connected_m{n,sd},smoothness))-min(smooth(t_connected_m{n,sd},smoothness))),  max(smooth(t_connected_m{n,sd},smoothness))+.2*(max(smooth(t_connected_m{n,sd},smoothness))-min(smooth(t_connected_m{n,sd},smoothness)))]);end
                    [hAx,hLine1,hLine2] = plotyy(tb, smooth(t_unconnected_m{n,sd}, 20), tb, smooth(t_connected_m{n,sd}-t_unconnected_m{n,sd}, 20));
                    diffcu =smooth(t_connected_m{n,sd}-t_unconnected_m{n,sd},20);
                    hAx(2).XLim = [-.2 0.9];
                    hLine1.LineWidth = 1;
                    hLine1.Color = [1 0 0];
                    try;hAx(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
                    hLine2.Color=[0 1 0];
                    hLine2.LineWidth = 1;
                    
                    
                    plot(tb, smooth(t_connected_m2{n,sd}, 20), 'b', 'LineWidth',1); ylim(ylim1)%try  ylim([min(smooth(t_connected_m{n,sd},smoothness))-(max(smooth(t_connected_m{n,sd},smoothness))-min(smooth(t_connected_m{n,sd},smoothness))),  max(smooth(t_connected_m{n,sd},smoothness))+.2*(max(smooth(t_connected_m{n,sd},smoothness))-min(smooth(t_connected_m{n,sd},smoothness)))]);end
                    [hAx,hLine1,hLine2] = plotyy(tb, smooth(t_unconnected_m2{n,sd}, 20), tb, smooth(t_connected_m2{n,sd}-t_unconnected_m2{n,sd}, 20));
                    diffcu =smooth(t_connected_m2{n,sd}-t_unconnected_m2{n,sd},20);
                    hAx(2).XLim = [-.2 0.9];
                    hLine1.LineWidth = 1;
                    hLine1.Color = [1 .5 .5];
                    try;hAx(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
                    hLine2.Color=[0 1 0];
                    hLine2.LineWidth = 1;
                    
                    
                    %             plot(tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, 20), 'g', 'LineWidth',2);
                end
                plot(tb, smooth(t_connected{n,sd}, 20),'b')
                [hAx2,hLine12,hLine22] = plotyy(tb, smooth(t_unconnected{n,sd}, 20),tb, smooth(t_connected{n,sd}-t_unconnected{n,sd}, 20)); xlim([-.2 0.9]);
                hLine12.Color = [1 0 0];
                hLine22.Color = [0 0 0];
                hAx2(2).XLim  = [-.2 0.9];
                axes(hAx2(2));line(hAx2(2).XLim,[0 0]);
                if ~plot_mask
                    diffcu =smooth(t_connected{n,sd}-t_unconnected{n,sd},20);
                end
                p=patch([1e-3*Poss_dur(sd) 1e-3*Poss_dur(sd)+.3 1e-3*Poss_dur(sd)+.3 1e-3*Poss_dur(sd)],[-10 -10 20 20],'r');
                set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
                
                try;hAx2(1).YLim = [-1 1.2*max(smooth(t_connected{n,sd}, 20))];hAx2(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
                %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
                
                
                
                
            end
            
                     
            
            AUC = zeros(1,numel(Poss_dur));
            OffResp = AUC;
            for sd = 1:numel(Poss_dur)
                t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
                diff = smooth(t_connected{n,sd}(t_relev)-t_unconnected{n,sd}(t_relev),20);
%                 diff(diff<0) = 0;
                AUC(sd) = sum(diff);
                tbl = find(tb>((1e-3*Poss_dur(sd))+0.01) & tb<((1e-3*Poss_dur(sd)) +0.03));
                blact = smooth(t_connected{n,sd}(tbl)-t_unconnected{n,sd}(tbl),20);
                pre_off = mean(blact);
                tpeak = find(tb>((1e-3*Poss_dur(sd))+0.03) & tb<((1e-3*Poss_dur(sd)) +0.09));
                pact = smooth(t_connected{n,sd}(tpeak)-t_unconnected{n,sd}(tpeak),20);
                peakvalue = max(pact);
                OffResp(sd) = peakvalue - pre_off;
                OffRatio(sd) = peakvalue/pre_off;
            end
            subplot(2,4,8);hold on; plot(Poss_dur, AUC, 'Color', [.1 .1 .1]);
            if plot_mask
                AUC_m = zeros(1,numel(Poss_dur));
                AUC_m2 = zeros(1,numel(Poss_dur));
                for sd = 1:numel(Poss_dur)
                    t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
                    diff_m =abs(smooth(t_connected_m{n,sd}(t_relev)-t_unconnected_m{n,sd}(t_relev), 20));
                    diff_m2 =abs( smooth(t_connected_m2{n,sd}(t_relev)-t_unconnected_m2{n,sd}(t_relev), 20));
%                     diff_m(diff_m<0)=0;
                    AUC_m(sd) = sum(diff_m);
                    AUC_m2(sd) = sum(diff_m2);
                end
                hold on; plot(Poss_dur, AUC_m, 'Color', [0 .9 0]);
                plot(Poss_dur, AUC_m2, 'Color', [.3 .9 .3]);
            end
            try
             flin = @(param,xval) param(1)+param(2)*xval;
             flin_m = @(param,xval) param(3)+param(2)*xval;
             flin_m2 = @(param,xval) param(4)+param(2)*xval;
             [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_dur, Poss_dur, Poss_dur}, {AUC,AUC_m,AUC_m2}, {flin, flin_m, flin_m2}, [20 0.2 20 20]);
             x2 = 0:10:Poss_dur(end);
             y2 = beta(1)+beta(2)*x2;
             y2_m =  beta(3)+beta(2)*x2;
             y2_m2 =  beta(4)+beta(2)*x2;
             plot(x2, y2, 'Color','k')
             plot(x2, y2_m, 'Color','g')
             plot(x2, y2_m2, 'Color',[.3 .9 .3])
             [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked 1frame', 'masked immediate', 'Location', 'NorthWest');
             LegendPos = get(legend_h,'Position');
            end
             EstimWorthNeuronal(n) = (beta(1)-beta(4))/beta(2);
             dim = [LegendPos(1) LegendPos(2)-LegendPos(4) .3 LegendPos(4)/2];
             str = ['Estimated worth: ' num2str(EstimWorthNeuronal(n),3) 'ms'];
             annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
             OffResp_chan(n) = mean(OffResp(end-2:end));
             OffRatio_cann(n) = mean(OffRatio(end-2:end));
             
             
            xlabel('stimulus duration')
            pause(0.02)
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-dpng',[SaveDir 'ConnUnconn_Chan' num2str(n) '_s' num2str(days(1)) '_' num2str(days(end)) '.png'])
            
        else
            figure(f2);clf
            
            hold on
            plot(tb,smooth(background{n,1},20),'k');
            plot(tb,smooth(ttree_t{n,1},20),'b');
            plot(tb,smooth(ttree_d{n,1},20),'r');
            plot(tb,smooth(dtree_t{n,1},20),'c');
            plot(tb,smooth(dtree_d{n,1},20),'m');
            legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'})
            plot(tb,smooth(background{n,2},20),'k', 'LineWidth',2);
            plot(tb,smooth(ttree_t{n,2},20),'b', 'LineWidth',2);
            plot(tb,smooth(ttree_d{n,2},20),'r', 'LineWidth',2);
            plot(tb,smooth(dtree_t{n,2},20),'c', 'LineWidth',2);
            plot(tb,smooth(dtree_d{n,2},20),'m', 'LineWidth',2); xlim([-.2 1.1]);
            title(num2str(n))
            pause(0.2)
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-dpng',[SaveDir 'AllcondNoDur_Chan' num2str(n) '_s' num2str(days(1)) '_' num2str(days(end)) '.png'])
            
            
            
            
            figure(f10);subplot(7,7,n)
            ylim1 = [-2 1.5];
            
            hold on; %title([ 'chan ' num2str(n)])
            
            if plot_mask
                plot(tb, smooth(t_connected_m{n}, 20), 'b', 'LineWidth',1); ylim(ylim1)%try  ylim([min(smooth(t_connected_m{n,sd},smoothness))-(max(smooth(t_connected_m{n,sd},smoothness))-min(smooth(t_connected_m{n,sd},smoothness))),  max(smooth(t_connected_m{n,sd},smoothness))+.2*(max(smooth(t_connected_m{n,sd},smoothness))-min(smooth(t_connected_m{n,sd},smoothness)))]);end
                [hAx,hLine1,hLine2] = plotyy(tb, smooth(t_unconnected_m{n}, 20), tb, smooth(t_connected_m{n}-t_unconnected_m{n}, 20));
                diffcu =smooth(t_connected_m{n,sd}-t_unconnected_m{n},20);
                hAx(2).XLim = [-.2 0.9];
                hLine1.LineWidth = 1;
                hLine1.Color = [1 0 0];
                try;hAx(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
                hLine2.Color=[0 1 0];
                hLine2.LineWidth = 1;
                
                %             plot(tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, 20), 'g', 'LineWidth',2);
            end
            plot(tb, smooth(t_connected{n}, 20),'b')
            [hAx2,hLine12,hLine22] = plotyy(tb, smooth(t_unconnected{n}, 20),tb, smooth(t_connected{n}-t_unconnected{n}, 20)); xlim([-.2 0.9]);
            hLine12.Color = [1 0 0];
            hLine22.Color = [0 0 0];
            hAx2(2).XLim  = [-.2 0.9];
            axes(hAx2(2));line(hAx2(2).XLim,[0 0]);
            if ~plot_mask
                diffcu =smooth(t_connected{n}-t_unconnected{n},20);
            end
            %                 p=patch([1e-3*Poss_dur(sd) 1e-3*Poss_dur(sd)+.3 1e-3*Poss_dur(sd)+.3 1e-3*Poss_dur(sd)],[-10 -10 20 20],'r');
            %                 set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
            hAx2(1).YLim = [-1 1.2*max(smooth(t_connected{n}, 20))];
            try;hAx2(2).YLim = [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))];end
            %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
            
            pause(0.02)
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-dpng',[SaveDir 'ConnUnconnNoDur_s' num2str(days(1)) '_' num2str(days(end)) '.png'])
            
            
            
        end
        
        
    end
    
    
    
    
end
%%
good_chann = find(SNR>SNR_threshold);
if plot_chann
    figure;plot(OffResp_chan(good_chann), EstimWorthNeuronal(good_chann), 'o')
    [R,P] = corrcoef(OffResp_chan(good_chann),EstimWorthNeuronal(good_chann));
    linear_fit = polyfit(OffResp_chan(good_chann),EstimWorthNeuronal(good_chann),1);
    x_fit = 0:0.1:(1.1*max(OffResp_chan(good_chann)));
    hold on; plot(x_fit, linear_fit(2)+linear_fit(1)*x_fit);
    xlabel('Amplitude of Offset Response')
    ylabel('Estimated Worth')
    title({'correlation between offset response and', 'neuronal worth at single channel level'}) 
    filename = ['CorrOffsetWorth_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
end
%%
if plot_stim_dur
    y_perf = perf(2,:)./sum(perf);
    y_perf_m = perf_m(2,:)./sum(perf_m);
    y_perf_m2 = perf_m2(2,:)./sum(perf_m2);
    
    figure('Position', [100 100 500 800], 'Color',[1 1 1]);  
    subplot(2,1,1);plot(Poss_dur, y_perf)
    if plot_mask
        hold on;  plot(Poss_dur, y_perf_m,'r')
         plot(Poss_dur, y_perf_m2,'Color', [1 .3 .3])
    end
    title('Behavioral performance on the iconic memory task');
    xlabel('stimulus duration (ms)')
    
    fsigm = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(2)-xval)*param(3)));
    fsigm_m = @(param,xval) 0.5+(param(4)-0.5)./(1+10.^((param(5)-xval)*param(6)));
    fsigm_m2 = @(param,xval) 0.5+(param(7)-0.5)./(1+10.^((param(8)-xval)*param(9)));

    [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_dur, Poss_dur, Poss_dur}, {y_perf,y_perf_m,y_perf_m}, {fsigm, fsigm_m, fsigm_m2}, [.8  150 0.005 .8  150 0.005 .8  150 0.005]);
    x1 = 0:10:Poss_dur(end);
    y1 = 0.5+(beta(1)-0.5)./(1+10.^((beta(2)-x1)*beta(3)));
    y1_m = 0.5+(beta(4)-0.5)./(1+10.^((beta(5)-x1)*beta(6)));
    y1_m2 = 0.5+(beta(7)-0.5)./(1+10.^((beta(8)-x1)*beta(9)));
    plot(x1, y1, 'Color','b')
    plot(x1, y1_m, 'Color','r')
    plot(x1, y1_m2, 'Color', [1 .3 .3])
    RSS2 = sum(r.^2);
    
 
    % plot with shared fit parameters
    subplot(2,1,2);plot(Poss_dur, y_perf)
    if plot_mask
        hold on;  plot(Poss_dur, y_perf_m,'r')
         plot(Poss_dur, y_perf_m2,'Color', [1 .3 .3])
    end
    
    fsigm = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(2)-xval)*param(3)));
    fsigm_m = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(4)-xval)*param(3)));
    fsigm_m2 = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(5)-xval)*param(3)));
    [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_dur, Poss_dur, Poss_dur}, {perf(2,:)./sum(perf),perf_m(2,:)./sum(perf_m), perf_m2(2,:)./sum(perf_m2)}, {fsigm, fsigm_m, fsigm_m2}, [.8  150 0.005 0.005 0.005]);
    RSS1 = sum(r.^2);
    F = ((RSS1-RSS2)/(9-5))/(RSS2/(numel(r)-9));
    p = fcdf(F, 9-5, (numel(r)-9));
    x1 = 0:10:Poss_dur(end);
    y1 = 0.5+(beta(1)-0.5)./(1+10.^((beta(2)-x1)*beta(3)));
    y1_m = 0.5+(beta(1)-0.5)./(1+10.^((beta(4)-x1)*beta(3)));
    y1_m2 = 0.5+(beta(1)-0.5)./(1+10.^((beta(5)-x1)*beta(3)));
    plot(x1, y1, 'Color','b')
    plot(x1, y1_m, 'Color','r')
    plot(x1, y1_m2, 'Color', [1 .3 .3])
    xlabel('stimulus duration (ms)')
    
    
    [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked 1f', 'masked imm','Location', 'SouthEast');
    LegendPos = get(legend_h,'Position');
    set(legend_h, 'Location', 'NorthWest')
    EstimeWorth1f = (beta(4)-beta(2));
    EstimWorthBehav = (beta(5)-beta(2));
    dim = [LegendPos(1)-LegendPos(3) LegendPos(2) .3 LegendPos(4)/2];
    dim1f = [LegendPos(1)-LegendPos(3) .05+LegendPos(2) .3 LegendPos(4)/2];
    str = ['Estimated worth: ' num2str(EstimWorthBehav,3) 'ms'];
    str1f = ['Estimated worth after 1frame: ' num2str(EstimeWorth1f,3) 'ms'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
    annotation('textbox',dim1f,'String',str1f,'FitBoxToText','on', 'LineStyle', 'none')
    dim_F = dim; dim_F(2) = dim_F(2)+2*LegendPos(4); dim_F(4) = LegendPos(4);
    if p<0.05
        str = {'The more complex model explains the data', ['significantly better  (p = ' num2str(p,3) ')']};
    else
        str = {'The more complex model does', 'not explain the data',  ['significantly better (p = ' num2str(p,3) ')']};
    end
    annotation('textbox',dim_F,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
    
    filename = ['performances_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    %%
    f1 = figure('Position', [100 50 1400 900]);
    if plot_mask; f2 = figure('Position', [100 50 1400 900]);end
    for sd = 1:min(7,numel(Poss_dur))
        t_connected_avg{sd} = 0*background{48,1,sd};n_t_connected_avg(sd)=0;
        t_unconnected_avg{sd} = 0*background{48,1,sd};n_t_unconnected_avg(sd)=0;
        for h = 1:2
            background_avg{h,sd} = 0*background{48,1,sd};n_background_avg(h,sd)=0;
            ttree_t_avg{h,sd} = 0*background{48,1,sd};n_ttree_t_avg(h,sd)=0;
            ttree_d_avg{h,sd} = 0*background{48,1,sd};n_ttree_d_avg(h,sd)=0;
            dtree_t_avg{h,sd} = 0*background{48,1,sd};n_dtree_t_avg(h,sd)=0;
            dtree_d_avg{h,sd} = 0*background{48,1,sd};n_dtree_d_avg(h,sd)=0;
        end
        if plot_mask
            t_connected_m_avg{sd} = 0*background{48,1,sd};n_t_connected_m_avg(sd)=0;
            t_unconnected_m_avg{sd} = 0*background{48,1,sd};n_t_unconnected_m_avg(sd)=0;
            for h = 2%1:
                background_m_avg{h,sd} = 0*background{48,1,sd};n_background_m_avg(h,sd)=0;
                ttree_t_m_avg{h,sd} = 0*background{48,1,sd};n_ttree_t_m_avg(h,sd)=0;
                ttree_d_m_avg{h,sd} = 0*background{48,1,sd};n_ttree_d_m_avg(h,sd)=0;
                dtree_t_m_avg{h,sd} = 0*background{48,1,sd};n_dtree_t_m_avg(h,sd)=0;
                dtree_d_m_avg{h,sd} = 0*background{48,1,sd};n_dtree_d_m_avg(h,sd)=0;
            end
            
            
            t_connected_m2_avg{sd} = 0*background{48,1,sd};n_t_connected_m2_avg(sd)=0;
            t_unconnected_m2_avg{sd} = 0*background{48,1,sd};n_t_unconnected_m2_avg(sd)=0;
            for h = 2%1:
                background_m2_avg{h,sd} = 0*background{48,1,sd};n_background_m2_avg(h,sd)=0;
                ttree_t_m2_avg{h,sd} = 0*background{48,1,sd};n_ttree_t_m2_avg(h,sd)=0;
                ttree_d_m2_avg{h,sd} = 0*background{48,1,sd};n_ttree_d_m2_avg(h,sd)=0;
                dtree_t_m2_avg{h,sd} = 0*background{48,1,sd};n_dtree_t_m2_avg(h,sd)=0;
                dtree_d_m2_avg{h,sd} = 0*background{48,1,sd};n_dtree_d_m2_avg(h,sd)=0;
            end
            
        end
        for n = good_chann
            
            try t_connected_avg{sd} = nansum([t_connected_avg{sd} ,t_connected{n,sd}],2);end; if ~isnan(sum(t_connected{n,sd}));n_t_connected_avg(sd) = n_t_connected_avg(sd)+1;end
            try t_unconnected_avg{sd} =  nansum([t_unconnected_avg{sd},t_unconnected{n,sd}],2);end;if ~isnan(sum(t_unconnected{n,sd}));n_t_unconnected_avg(sd) = n_t_unconnected_avg(sd)+1;end
            for h = 1:2
                try background_avg{h,sd} =  nansum([background_avg{h,sd} ,background{n,h,sd}],2);end;if ~isnan(sum(background{n,h,sd}));n_background_avg(h,sd) = n_background_avg(h,sd)+1;end
                try ttree_t_avg{h,sd} =  nansum([ttree_t_avg{h,sd} ,ttree_t{n,h,sd}],2);end ;if ~isnan(sum(ttree_t{n,h,sd}));n_ttree_t_avg(h,sd) = n_ttree_t_avg(h,sd)+1;end
                try ttree_d_avg{h,sd} =  nansum([ttree_d_avg{h,sd} ,ttree_d{n,h,sd} ],2);end;if ~isnan(sum(ttree_d{n,h,sd}));n_ttree_d_avg(h,sd) = n_ttree_d_avg(h,sd)+1;end
                
                try dtree_t_avg{h,sd} =  nansum([dtree_t_avg{h,sd},dtree_t{n,h,sd}],2);end;if ~isnan(sum(dtree_t{n,h,sd}));n_dtree_t_avg(h,sd) = n_dtree_t_avg(h,sd)+1;end
                try dtree_d_avg{h,sd} =  nansum([dtree_d_avg{h,sd},dtree_d{n,h,sd}],2);end;if ~isnan(sum(dtree_d{n,h,sd}));n_dtree_d_avg(h,sd) = n_dtree_d_avg(h,sd)+1;end
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
                
                try t_connected_m2_avg{sd} =  nansum([t_connected_m2_avg{sd} ,t_connected_m2{n,sd}],2);end;if ~isnan(sum(t_connected_m2{n,sd}));n_t_connected_m2_avg(sd) = n_t_connected_m2_avg(sd)+1;end
                try t_unconnected_m2_avg{sd} =  nansum([t_unconnected_m2_avg{sd} ,t_unconnected_m2{n,sd}],2);end;if ~isnan(sum(t_unconnected_m2{n,sd}));n_t_unconnected_m2_avg(sd) = n_t_unconnected_m2_avg(sd)+1;end
                for h = 1:2
                    try background_m2_avg{h,sd} =  nansum([background_m2_avg{h,sd} ,background_m2{n,h,sd}],2);end;if ~isnan(sum(background_m2{n,h,sd}));n_background_m2_avg(h,sd) = n_background_m2_avg(h,sd)+1;end
                    try  ttree_t_m2_avg{h,sd} =  nansum([ttree_t_m2_avg{h,sd} ,ttree_t_m2{n,h,sd}],2);end;if ~isnan(sum(ttree_t_m2{n,h,sd}));n_ttree_t_m2_avg(h,sd) = n_ttree_t_m2_avg(h,sd)+1;end
                    try ttree_d_m2_avg{h,sd} =  nansum([ttree_d_m2_avg{h,sd} ,ttree_d_m2{n,h,sd} ],2);end;if ~isnan(sum(ttree_d_m2{n,h,sd}));n_ttree_d_m2_avg(h,sd) = n_ttree_d_m2_avg(h,sd)+1;end
                    
                    try dtree_t_m2_avg{h,sd} =  nansum([dtree_t_m2_avg{h,sd},dtree_t_m2{n,h,sd}],2);end;if ~isnan(sum(dtree_t_m2{n,h,sd}));n_dtree_t_m2_avg(h,sd) = n_dtree_t_m2_avg(h,sd)+1;end
                    try dtree_d_m2_avg{h,sd} =  nansum([dtree_d_m2_avg{h,sd},dtree_d_m2{n,h,sd}],2);end;if ~isnan(sum(dtree_d_m2{n,h,sd}));n_dtree_d_m2_avg(h,sd) = n_dtree_d_m2_avg(h,sd)+1;end
                end
                
                
            end
        end
        
        t_connected_avg{sd} = t_connected_avg{sd}/n_t_connected_avg(sd);
        t_unconnected_avg{sd} = t_unconnected_avg{sd}/n_t_unconnected_avg(sd);
        for h = 1:2
            background_avg{h,sd} = background_avg{h,sd}/n_background_avg(h,sd);
            ttree_t_avg{h,sd} = ttree_t_avg{h,sd}/n_ttree_t_avg(h,sd);
            ttree_d_avg{h,sd} = ttree_d_avg{h,sd}/n_ttree_d_avg(h,sd);
            dtree_t_avg{h,sd} = dtree_t_avg{h,sd}/n_dtree_t_avg(h,sd);
            dtree_d_avg{h,sd} = dtree_d_avg{h,sd}/n_dtree_d_avg(h,sd);
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
            
            
            t_connected_m2_avg{sd} = t_connected_m2_avg{sd}/n_t_connected_m2_avg(sd);
            t_unconnected_m2_avg{sd} = t_unconnected_m2_avg{sd}/n_t_unconnected_m2_avg(sd);
            for h = 1:2
                background_m2_avg{h,sd} = background_m2_avg{h,sd}/n_background_m2_avg(h,sd);
                ttree_t_m2_avg{h,sd} = ttree_t_m2_avg{h,sd}/n_ttree_t_m2_avg(h,sd);
                ttree_d_m2_avg{h,sd} = ttree_d_m2_avg{h,sd}/n_ttree_d_m2_avg(h,sd);
                dtree_t_m2_avg{h,sd} = dtree_t_m2_avg{h,sd}/n_dtree_t_m2_avg(h,sd);
                dtree_d_m2_avg{h,sd} = dtree_d_m2_avg{h,sd}/n_dtree_d_m2_avg(h,sd);
            end
        end
        
        
        figure(f1);subplot(2,4,sd);hold on; title([num2str(Poss_dur(sd)), ' ms'])%
        %         plot(tb,smooth(background_avg{1,sd},20),'k');
        %         plot(tb,smooth(ttree_t_avg{1,sd},20),'b');
        %         plot(tb,smooth(ttree_d_avg{1,sd},20),'r');
        %         plot(tb,smooth(dtree_t_avg{1,sd},20),'c');
        %         plot(tb,smooth(dtree_d_avg{1,sd},20),'m');
        plot(tb,smooth(background_avg{2,sd},20),'k', 'LineWidth',2);
        plot(tb,smooth(ttree_t_avg{2,sd},20),'b', 'LineWidth',2);
        plot(tb,smooth(ttree_d_avg{2,sd},20),'r', 'LineWidth',2);
        plot(tb,smooth(dtree_t_avg{2,sd},20),'c', 'LineWidth',2);
        plot(tb,smooth(dtree_d_avg{2,sd},20),'m', 'LineWidth',2); xlim([-.2 0.9]);
        if sd == 1; legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'});end
        
        if plot_mask
            figure(f2);subplot(2,4,sd);hold on; title([num2str(Poss_dur(sd)), ' ms, masked at 20ms'])%
            %             plot(tb,smooth(background_m_avg{1,sd},20),'k');
            %             plot(tb,smooth(ttree_t_m_avg{1,sd},20),'b');
            %             plot(tb,smooth(ttree_d_m_avg{1,sd},20),'r');
            %             plot(tb,smooth(dtree_t_m_avg{1,sd},20),'c');
            %             plot(tb,smooth(dtree_d_m_avg{1,sd},20),'m');
            plot(tb,smooth(background_m_avg{2,sd},20),'k', 'LineWidth',2);
            plot(tb,smooth(ttree_t_m_avg{2,sd},20),'b', 'LineWidth',2);
            plot(tb,smooth(ttree_d_m_avg{2,sd},20),'r', 'LineWidth',2);
            plot(tb,smooth(dtree_t_m_avg{2,sd},20),'c', 'LineWidth',2);
            plot(tb,smooth(dtree_d_m_avg{2,sd},20),'m', 'LineWidth',2); xlim([-.2 0.9]);
            
            
                        
            plot(tb,smooth(background_m2_avg{2,sd},20),'k', 'LineWidth',3);
            plot(tb,smooth(ttree_t_m2_avg{2,sd},20),'b', 'LineWidth',3);
            plot(tb,smooth(ttree_d_m2_avg{2,sd},20),'r', 'LineWidth',3);
            plot(tb,smooth(dtree_t_m2_avg{2,sd},20),'c', 'LineWidth',3);
            plot(tb,smooth(dtree_d_m2_avg{2,sd},20),'m', 'LineWidth',3); xlim([-.2 0.9]);
            if sd == 1; legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'});end
            
            
        end
    end
    
    
    
    figure(f1)
    filename = 'NeurActiv_withStimDur_noMask';
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    if plot_mask
        figure(f2)
        filename = ['NeurActiv_withStimDur_Mask_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
        set(gcf, 'PaperPositionMode', 'auto')
        print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    end
    %%
    figure('Position', [100 50 1400 900])
    
    for sd = 1:min(numel(Poss_dur),7)
        subplot(2,4,sd);hold on; title(['Stim duration: ' num2str(Poss_dur(sd)), ' ms']); xlabel('time(s)'); ylabel('normalized activity')
        ylim1 = [-2 2];
        if plot_mask
            a = plot(tb, smooth(t_connected_m_avg{sd}, 20), 'b', 'LineWidth',1); ylim(ylim1)
            [hAx,hLine1,hLine2] = plotyy(tb, smooth(t_unconnected_m_avg{sd}, 20), tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, 20));
            hAx(2).XLim = [-.2 0.9];
            hLine1.LineWidth = 1;
            hLine1.Color = [1 0 0];
            hAx(1).YLim = ylim1;
            hAx(2).YLim = [-.6 3];
            hLine2.Color=[0 1 0];
            hLine2.LineWidth = 1;
            
            a2 = plot(tb, smooth(t_connected_m2_avg{sd}, 20), 'b', 'LineWidth',1); ylim(ylim1)
            [hAxm2,hLine1m2,hLine2m2] = plotyy(tb, smooth(t_unconnected_m2_avg{sd}, 20), tb, smooth(t_connected_m2_avg{sd}-t_unconnected_m2_avg{sd}, 20));
            hAxm2(2).XLim = [-.2 0.9];
            hLine1m2.LineWidth = 1;
            hLine1m2.Color = [1 .3 .3];
            hAxm2(1).YLim = ylim1;
            hAxm2(2).YLim = [-.6 3];
            hLine2m2.Color=[.3 1 .3];
            hLine2m2.LineWidth = 1;
            
            
            if sd==2
                legend([a, hLine1, hLine2,a2, hLine1m2, hLine2m2,], {'masked connected 1f', 'masked unconnected 1f', 'masked 1f conn-unconn' 'masked connected 0f', 'masked unconnected 0f', 'masked 0f conn-unconn'})
            end
            %             plot(tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, 20), 'g', 'LineWidth',2);
        end
        
        b = plot(tb, smooth(t_connected_avg{sd}, 20),'Color',[0 0 .7]);  ylim(ylim1)
        [hAx2,hLine12,hLine22] = plotyy(tb, smooth(t_unconnected_avg{sd}, 20),tb, smooth(t_connected_avg{sd}-t_unconnected_avg{sd}, 20)); xlim([-.2 0.9]);
        hLine12.Color = [.7 0 0];
        hLine22.Color=[0 0 0];
        hAx2(2).XLim = [-.2 0.9];
        hAx2(2).YLim = [-.6 3];
        hAx2(1).YLim = ylim1;
        axes(hAx2(2));line(hAx2(2).XLim,[0 0])
        p=patch([0.01 1e-3*max(Poss_dur)+integ_dur 1e-3*max(Poss_dur)+integ_dur 0.01],[-10 -10 20 20],'r');
        set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
        if sd==1
            legend([b, hLine12, hLine22], {'no mask connected', 'no mask unconnected', 'conn - unconn'})
        end
        %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
    end
    
    subplot(2,4,8);
    %     integ_dur = 0.2;
    smoothness = 30;
    AUC = zeros(1,min(numel(Poss_dur),7));
%     AUC2 = zeros(1,5);
    for sd = 1:min(numel(Poss_dur),7)
        t_relev = find(tb>(0.01) & tb<((1e-3*max(Poss_dur))+integ_dur));
        diff_m =abs(smooth(t_connected_avg{sd}, smoothness)-smooth(t_unconnected_avg{sd}, smoothness));
        diff_m = diff_m(t_relev);%diff_m(diff_m<0)=0;
        %                     i = find(diff_m<0, 1);diff_m(i:end)=0;
        AUC(sd) = sum(diff_m);
        %             AUC_m(sd) = sum(smooth(t_connected_m_avg{sd}(t_relev)-t_unconnected_m_avg{sd}(t_relev), smoothness));
    end
    hold on; plot(Poss_dur(1:min(numel(Poss_dur),7)), AUC, 'k')
    xlabel('Stimulus duration(ms)')
    
    
    
    AUC_m = zeros(1,min(numel(Poss_dur),7));
    AUC_m2 = zeros(1,min(numel(Poss_dur),7));
    for sd = 1:min(numel(Poss_dur),7)
        t_relev = find(tb>(0.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
        diff_m =abs(smooth(t_connected_m_avg{sd}, smoothness)-smooth(t_unconnected_m_avg{sd}, smoothness));
        diff_m = diff_m(t_relev);%diff_m(diff_m<0)=0;
        %                     i = find(diff_m<0, 1);diff_m(i:end)=0;
        AUC_m(sd) = sum(diff_m);
        
        
        diff_m2 =abs(smooth(t_connected_m2_avg{sd}, smoothness)-smooth(t_unconnected_m2_avg{sd}, smoothness));
        diff_m2 = diff_m2(t_relev);%diff_m(diff_m<0)=0;
        AUC_m2(sd) = sum(diff_m2);

        %             AUC_m(sd) = sum(smooth(t_connected_m_avg{sd}(t_relev)-t_unconnected_m_avg{sd}(t_relev), smoothness));
    end
    hold on; plot(Poss_dur(1:min(numel(Poss_dur),7)), AUC_m, 'g'); plot(Poss_dur(1:min(numel(Poss_dur),7)), AUC_m2,'Color', [.3 1 .3]);
    title({'Integral of the difference in activity', ['between stim onset and ' num2str(integ_dur*1000) 'ms after the end of the stimulus']})
    legend({'without mask', 'mask 1frame',  'imm mask'}, 'Location', 'NorthWest')
    C = mean(AUC) - mean(AUC_m);
    C2 = mean(AUC) - mean(AUC_m2);
    
    filename = ['Diff_ConnectedUnconn_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    
    %     figure('Position', [100 50 1400 900])
    %     integ_dur = 0.3;
    %     for sd = 1:min(numel(Poss_dur),7)
    %         subplot(2,4,sd);hold on; title([num2str(Poss_dur(sd)), ' ms'])
    %         ylim1 = [-1.5 1.5];
    %         if plot_mask
    %             plot(tb, smooth(t_connected_m_avg{sd}, 20), 'b', 'LineWidth',1); ylim(ylim1)
    %             [hAx,hLine1,hLine2] = plotyy(tb, smooth(t_unconnected_m_avg{sd}, 20), tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, 20));
    %             hAx(2).XLim = [-.2 0.9];
    %             hLine1.LineWidth = 1;
    %             hLine1.Color = [1 0 0];
    %             hAx(1).YLim = ylim1;
    %             hAx(2).YLim = [-.6 4];
    %             hLine2.Color=[0 1 0];
    %             hLine2.LineWidth = 1;
    %             %             plot(tb, smooth(t_connected_m_avg{sd}-t_unconnected_m_avg{sd}, 20), 'g', 'LineWidth',2);
    %         end
    %         plot(tb, smooth(t_connected_avg{sd}, 20),'b');  ylim(ylim1)
    %         [hAx2,hLine12,hLine22] = plotyy(tb, smooth(t_unconnected_avg{sd}, 20),tb, smooth(t_connected_avg{sd}-t_unconnected_avg{sd}, 20)); xlim([-.2 0.9]);
    %         hLine12.Color = [1 0 0];
    %         hLine22.Color=[0 0 0];
    %         hAx2(2).XLim = [-.2 0.9];
    %         hAx2(2).YLim = [-.6 4];
    %         hAx2(1).YLim = ylim1;
    %         axes(hAx2(2));line(hAx2(2).XLim,[0 0])
    %         p=patch([0.01 1e-3*Poss_dur(sd)+integ_dur 1e-3*Poss_dur(sd)+integ_dur 0.01],[-10 -10 20 20],'r');
    %         set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
    %         %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
    %     end
    %
    %     subplot(2,4,8);
    %     %     integ_dur = 0.2;
    %     smoothness = 30;
    %     AUC = zeros(1,5);
    %     for sd = 1:5
    %         t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
    %         diff_m =smooth(t_connected_avg{sd}, smoothness)-smooth(t_unconnected_avg{sd}, smoothness);
    %         diff_m = diff_m(t_relev);%diff_m(diff_m<0)=0;
    %         %                     i = find(diff_m<0, 1);diff_m(i:end)=0;
    %         AUC(sd) = sum(diff_m);
    %         %             AUC_m(sd) = sum(smooth(t_connected_m_avg{sd}(t_relev)-t_unconnected_m_avg{sd}(t_relev), smoothness));
    %     end
    %     hold on; plot(Poss_dur(1:5), AUC, 'k')
    %
    %
    %
    %
    %     AUC_m = zeros(1,5);
    %     for sd = 1:5
    %         t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
    %         diff_m =smooth(t_connected_m_avg{sd}, smoothness)-smooth(t_unconnected_m_avg{sd}, smoothness);
    %         diff_m = diff_m(t_relev);%diff_m(diff_m<0)=0;
    %         %                     i = find(diff_m<0, 1);diff_m(i:end)=0;
    %         AUC_m(sd) = sum(diff_m);
    %         %             AUC_m(sd) = sum(smooth(t_connected_m_avg{sd}(t_relev)-t_unconnected_m_avg{sd}(t_relev), smoothness));
    %     end
    %     hold on; plot(Poss_dur(1:5), AUC_m, 'g')
    %     title(['Integral of the difference in activity until ' num2str(integ_dur*1000) 'ms after the end of the stimulus'])
    %     legend({'without mask', 'with mask'})
    %     C = mean(AUC) - mean(AUC_m);
    %
    %     filename = ['Diff_ConnectedUnconn_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    %     set(gcf, 'PaperPositionMode', 'auto')
    %     print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    %
    
    %% Compute sum of difference in the "masked condition" for the whole stim duration for all stim durations
    Stim_info = zeros(1,numel(Poss_dur));
    for sd = 1:min(numel(Poss_dur),7)
        t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
        diff_m =abs(smooth(t_connected_avg{sd}, smoothness)-smooth(t_unconnected_avg{sd}, smoothness));
        diff_m = diff_m(t_relev);
        Stim_info(sd) = sum(diff_m);
    end
    figure;plot(Poss_dur, Stim_info, 'Color', [.1 .1 .1]);
   
    Stim_info_m = zeros(1,numel(Poss_dur));
    for sd = 1:min(numel(Poss_dur),7)
        t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
        diff_m =abs(smooth(t_connected_m_avg{sd}, smoothness)-smooth(t_unconnected_m_avg{sd}, smoothness));
        diff_m = diff_m(t_relev);
        Stim_info_m(sd) = sum(diff_m);
    end
    
    Stim_info_m2 = zeros(1,numel(Poss_dur));
    for sd = 1:min(numel(Poss_dur),7)
        t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
        diff_m2 =abs(smooth(t_connected_m2_avg{sd}, smoothness)-smooth(t_unconnected_m2_avg{sd}, smoothness));
        diff_m2 = diff_m2(t_relev);
        Stim_info_m2(sd) = sum(diff_m2);
    end
    
    hold on;plot(Poss_dur, Stim_info_m, 'Color', [0 .9 0]);plot(Poss_dur, Stim_info_m2, 'Color', [.3 .9 .3]);
    
    flin = @(param,xval) param(1)+param(2)*xval;
    flin_m = @(param,xval) param(3)+param(2)*xval;
    flin_m2 = @(param,xval) param(4)+param(2)*xval;
    [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_dur, Poss_dur, Poss_dur}, {Stim_info,Stim_info_m,Stim_info_m2}, {flin, flin_m, flin_m2}, [20 0.2 20 20 ]);
    x2 = 0:10:Poss_dur(end);
    y2 = beta(1)+beta(2)*x2;
    y2_m =  beta(3)+beta(2)*x2;
    y2_m2 =  beta(4)+beta(2)*x2;
    plot(x2, y2, 'Color','k')
    plot(x2, y2_m, 'Color','g')
    plot(x2, y2_m2, 'Color', [.3 .9 .3]);
    [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked 1f', 'masked 0f', 'Location', 'NorthWest');
    LegendPos = get(legend_h,'Position');
    EstimWorthNeuronal = (beta(1)-beta(4))/beta(2);
    EstimWorthNeuronal1f = (beta(1)-beta(3))/beta(2);
    dim = [LegendPos(1) LegendPos(2)-LegendPos(4) .3 LegendPos(4)/2];
    dim1f = [LegendPos(1) LegendPos(2)-2*LegendPos(4) .3 LegendPos(4)/2];
    str = ['Estimated worth: ' num2str(EstimWorthNeuronal,3) 'ms'];
    str1f = ['Estimated worth after 1 frame: ' num2str(EstimWorthNeuronal1f,3) 'ms'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
    annotation('textbox',dim1f,'String',str1f,'FitBoxToText','on', 'LineStyle', 'none')
    title('Information about the connections in function of stimulus duration')
    xlabel('stimulus duration')
    
    filename = ['Info_stimDur' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    %%
    % compute integral of diff between connected and unconn curves over the
    % interval 50-400
    %     t_relev = find(tb>.01 & tb<0.5);
    if 0
    AUC = zeros(1,numel(Poss_dur));
    for sd = 1:min(numel(Poss_dur),7)
        t_relev = find(tb>(Poss_dur(sd)/1e3) & tb<(Poss_dur(sd)/1e3 +.3));
        AUC(sd) = sum(abs(smooth(t_connected_avg{sd}(t_relev)-t_unconnected_avg{sd}(t_relev),20)));
    end
    figure; plot(Poss_dur, AUC)
    if plot_mask
        AUC_m = zeros(1,numel(Poss_dur));
        for sd = 1:min(numel(Poss_dur),7)
            t_relev = find(tb>(Poss_dur(sd)/1e3) & tb<(Poss_dur(sd)/1e3 +.3));
            AUC_m(sd) = sum(abs(smooth(t_connected_m_avg{sd}(t_relev)-t_unconnected_m_avg{sd}(t_relev), 20)));
        end
        hold on; plot(Poss_dur, AUC_m, 'r')
    end
    legend({'no mask', 'masked'})
    xlabel('stimulus duration')
    title('Area Under the difference curve from stim end to stim end+300ms')
    filename = 'Diff_ConnectedUnconn_AUC';
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    
    sd = 3;
    a =100:10:500;
    cum_int = zeros(1, numel(a));
    for t = 1:numel(a)
        t_relev = find(tb>.01 & tb<a(t)/1000);
        cum_int(t) = sum(t_connected_avg{sd}(t_relev)-t_unconnected_avg{sd}(t_relev));
    end
    figure; plot(a, cum_int)
    %% 100ms of information :AUC
    t_relev = find(tb>0.1 & tb<.2);
    value = .5*(sum(t_connected_avg{5}(t_relev)-t_unconnected_avg{5}(t_relev))+ sum(t_connected_avg{4}(t_relev)-t_unconnected_avg{4}(t_relev))) ;%- sum(t_connected_avg{3}(t_relev)-t_unconnected_avg{3}(t_relev));
    end
    %%
    
else
    figure;
    for h = 1:2
        background_avg{h} = 0*background{48,h} ;n_background_avg(h)=0;
        ttree_t_avg{h} = 0*background{48,h};n_ttree_t_avg(h)=0;
        ttree_d_avg{h} = 0*background{48,h};n_ttree_d_avg(h)=0;
        dtree_t_avg{h} = 0*background{48,h};n_dtree_t_avg(h)=0;
        dtree_d_avg{h} = 0*background{48,h};n_dtree_d_avg(h)=0;
    end
    if plot_mask
        for h = 1:2
            background_m_avg{h} = 0*background{48,h} ;n_background_m_avg(h)=0;
            ttree_t_m_avg{h} = 0*background{48,h};n_ttree_t_m_avg(h)=0;
            ttree_d_m_avg{h} = 0*background{48,h};n_ttree_d_m_avg(h)=0;
            dtree_t_m_avg{h} = 0*background{48,h};n_dtree_t_m_avg(h)=0;
            dtree_d_m_avg{h} = 0*background{48,h};n_dtree_d_m_avg(h)=0;
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
    plot(tb,smooth(background_avg{1},20),'k');
    plot(tb,smooth(ttree_t_avg{1},20),'b');
    plot(tb,smooth(ttree_d_avg{1},20),'r');
    plot(tb,smooth(dtree_t_avg{1},20),'c');
    plot(tb,smooth(dtree_d_avg{1},20),'m');
    legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'})
    plot(tb,smooth(background_avg{2},20),'k', 'LineWidth',2);
    plot(tb,smooth(ttree_t_avg{2},20),'b', 'LineWidth',2);
    plot(tb,smooth(ttree_d_avg{2},20),'r', 'LineWidth',2);
    plot(tb,smooth(dtree_t_avg{2},20),'c', 'LineWidth',2);
    plot(tb,smooth(dtree_d_avg{2},20),'m', 'LineWidth',2); xlim([-.2 1.1]);
    if plot_mask
        subplot(1,2,2); hold on
        plot(tb,smooth(background_m_avg{1},20),'k');
        plot(tb,smooth(ttree_t_m_avg{1},20),'b');
        plot(tb,smooth(ttree_d_m_avg{1},20),'r');
        plot(tb,smooth(dtree_t_m_avg{1},20),'c');
        plot(tb,smooth(dtree_d_m_avg{1},20),'m');
        legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'})
        plot(tb,smooth(background_m_avg{2},20),'k', 'LineWidth',2);
        plot(tb,smooth(ttree_t_m_avg{2},20),'b', 'LineWidth',2);
        plot(tb,smooth(ttree_d_m_avg{2},20),'r', 'LineWidth',2);
        plot(tb,smooth(dtree_t_m_avg{2},20),'c', 'LineWidth',2);
        plot(tb,smooth(dtree_d_m_avg{2},20),'m', 'LineWidth',2); xlim([-.2 1.1]);
        
    end
end


% figure,bar(SNR)
% figure,bar(peak)
%
% figure; for t = 1:155; plot(tb, e(:,t)); pause(0.1);end
