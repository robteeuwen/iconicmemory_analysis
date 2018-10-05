% analysis function for the worth/decay of iconic memory estimate 
% This script can be used for the 'control' experiment in Monty with both
% variable stimulus duration and variable mask onset. Sessions 33 to 47
% written by Catherine Wacongne (catherine.waco@gmail.com) 


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
Sessions = 37:47;
for n = 1:48
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
                mask = mask.*(LOG.Masked(trials)) +1;
                mask_delay = LOG.Mask_delay(trials);
                mask_delay(mask==1) = 500;
                Poss_delays = unique(mask_delay);
                
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
                mask =  [mask  mask_d.*(LOG.Masked(trials))+1];
                mask_delay = [mask_delay LOG.Mask_delay(trials)];
                mask_delay(mask==1) = 500;
                Poss_delays = [Poss_delays unique(mask_delay)];
                
            end
            
        end
    end
    Poss_dur = unique(Poss_dur);
    Poss_delays = unique(Poss_delays);
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
        for sd = 1:numel(Poss_dur)
            for md = [1:numel(Poss_delays)]
                % average per condition
                ind = find(isdrawn4==0 & Hit==h & stim_dur == Poss_dur(sd) & mask_delay==Poss_delays(md) );%& stim_dur == Poss_dur(sd) & mask == m
                background{n,h,sd,md}  = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(targ==10 & Hit==h & stim_dur == Poss_dur(sd) & mask_delay==Poss_delays(md) );
                ttree_t{n,h,sd,md} = nanmean(AllTrials(:,ind),2);
                
                
                
                ind = find(targ==9 & Hit==h & stim_dur == Poss_dur(sd) & mask_delay==Poss_delays(md) );
                ttree_d{n,h,sd,md} = nanmean(AllTrials(:,ind),2);
                
                ind = find(connected==1 & targ~=9 & Hit==h & stim_dur == Poss_dur(sd)& mask_delay==Poss_delays(md)  );
                dtree_t{n,h,sd,md} = nanmean(AllTrials(:,ind),2);
                
                
                ind = find(connected==2 & targ~=10 & Hit==h & stim_dur == Poss_dur(sd)& mask_delay==Poss_delays(md)  );
                dtree_d{n,h,sd,md} = nanmean(AllTrials(:,ind),2);
                
                ind = find( Hit==h & stim_dur == Poss_dur(sd) & mask_delay==Poss_delays(md) );
                perf(h,sd,md) = numel(ind);
                
                
                ind = find(connected==1 & stim_dur == Poss_dur(sd) & mask_delay==Poss_delays(md) );% & Hit==2);%
                t_connected{n,sd,md} =  nanmean(AllTrials(:,ind),2);
                ind = find(connected==2 & stim_dur == Poss_dur(sd) & mask_delay==Poss_delays(md) );%& Hit==2);%
                t_unconnected{n,sd,md} = nanmean(AllTrials(:,ind),2);
                
                
            end
        end
        
    end
    
end

%% plot the conditions

close all
if plot_chann
    f2 =figure('Position', [100 50 1400 900]);
    if ~plot_stim_dur
        f10 = figure('Position', [100 50 1400 900]);
    end
    for n = 1:48 
       
            figure(f2);clf
            for sd = 1:min(numel(Poss_dur),7)
                ylim1 = [-2 1.5];
                
                
                for md = [1:numel(Poss_delays)]
                    figure(f2)
                    subplot(2,4,sd);hold on; title([num2str(Poss_dur(sd)), ' ms, chan ' num2str(n)])
                   
                    yyaxis left
                    plot(tb, smooth(t_connected{n,sd,md}, 20),'Color', [1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays)) 1], 'LineStyle', '-', 'Marker','none' );hold on
                    plot(tb, smooth(t_unconnected{n,sd,md}, 20),'Color', [1 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))] , 'LineStyle', '-', 'Marker','none');
                    yyaxis right
                    plot(tb, smooth(t_connected{n,sd,md}-t_unconnected{n,sd,md}, 20),'Color',[1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))], 'LineStyle', '-', 'Marker','none'); xlim([-.2 0.9]);hold on
                    diffcu =smooth(t_connected{n,sd,md}-t_unconnected{n,sd,md},20);
%                   
                    
                    
                end
                p=patch([0.01 1e-3*Poss_dur(sd)+integ_dur 1e-3*Poss_dur(sd)+integ_dur 0.01],[-10 -10 20 20],'r');
                set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
                try
                yyaxis left 
                set(gca, 'YLim', [-1 1.2*max(smooth(t_connected{n,sd,md}, 20))])
                yyaxis right
                set(gca, 'YLim', [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))])
%                 keyboard
                end
                
            end
            
            
            
            AUC = zeros(numel(Poss_delays),numel(Poss_dur));
            OffResp = AUC;
            for md = [1:numel(Poss_delays)]
                for sd = 1:numel(Poss_dur)
                    t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
                    diff = abs(smooth(t_connected{n,sd,md}(t_relev)-t_unconnected{n,sd,md}(t_relev),20));
                    
                    AUC(md,sd) = sum(diff);
                    tbl = find(tb>((1e-3*Poss_dur(sd))+0.01) & tb<((1e-3*Poss_dur(sd)) +0.03));
                    blact = smooth(t_connected{n,sd,md}(tbl)-t_unconnected{n,sd,md}(tbl),20);
                    pre_off = mean(blact);
                    tpeak = find(tb>((1e-3*Poss_dur(sd))+0.03) & tb<((1e-3*Poss_dur(sd)) +0.09));
                    pact = smooth(t_connected{n,sd,md}(tpeak)-t_unconnected{n,sd,md}(tpeak),20);
                    peakvalue = max(pact);
                    OffResp(md,sd) = peakvalue - pre_off;
                    OffRatio(md,sd) = peakvalue/pre_off;
                   

                end
                OffResp_chan(md,n) = mean(OffResp(md,end-2:end));
                OffRatio_cann(md,n) = mean(OffRatio(md,end-2:end));
                 subplot(2,4,8);hold on; plot(Poss_dur, AUC(md,:), 'Color', [1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))]);
            end
            
            try
                clear Xreg flin Yreg
                gestimates = 0.2;
                for md = [1:numel(Poss_delays)-1]
                    Xreg{md} = Poss_dur;
                    str = ['@(param,xval)  param(' num2str(md+1) ')+param(1)*xval'];  
                    flin_m = str2func(str);
                    flin{md} = flin_m;
                    Yreg{md} = AUC(md,:);
                    gestimates = [gestimates 20];
                end
%                 @(param,xval) param(1)+param(2)*xval;
                
                [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(Xreg, Yreg, flin, gestimates);
                x2 = 0:10:Poss_dur(end);
                for md = 1:numel(Poss_delays)-1
                    y2 = beta(md+1)+beta(1)*x2;
                    plot(x2, y2, 'Color',[1-(md/numel(Poss_delays)) 1 1-(md/numel(Poss_delays))])
                    EstimWorthNeuronal(md,n) = (beta(end)-beta(md+1))/beta(1);
                end
                [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked', 'Location', 'NorthWest');           
                LegendPos = get(legend_h,'Position');
                set(legend_h,'visible','off')
            
            
            
                dim = [LegendPos(1) LegendPos(2)-LegendPos(4) .3 LegendPos(4)/2];
                str = ['Estimated worth: ' num2str(EstimWorthNeuronal(1,n),3) 'ms'];
                annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
                
            catch 
                warning('no neuronal worth was computed')
%                 keyboard
            end
            
            xlabel('stimulus duration')
            pause(0.02)
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-dpng',[SaveDir 'ConnUnconn_Chan' num2str(n) '_s' num2str(days(1)) '_' num2str(days(end)) '.png'])
            
       
        
    end
    
    
    
    
end
%%

%good_chann = [74:76 78:83 86 87 89 90:48];%74:48;%[75 81 82 83 85:89 91:48];%[74 75 78 83 86 90:48];
good_chann = find(SNR>SNR_threshold);
if plot_chann
    try
    figure;plot(OffResp_chan(end,good_chann), EstimWorthNeuronal(good_chann), 'o')
    [R,P] = corrcoef(OffResp_chan(end,good_chann),EstimWorthNeuronal(good_chann));
    linear_fit = polyfit(OffResp_chan(end,good_chann),EstimWorthNeuronal(good_chann),1);
    x_fit = 0:0.1:(1.1*max(OffResp_chan(end,good_chann)));
    hold on; plot(x_fit, linear_fit(2)+linear_fit(1)*x_fit);
    xlabel('Amplitude of Offset Response')
    ylabel('Estimated Worth')
    title({'correlation between offset response and', 'neuronal worth at single channel level'})
    filename = ['CorrOffsetWorth_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    end
end
%%
if plot_stim_dur
    y_perf = squeeze(perf(2,:,:)./sum(perf));
%     y_perf_m = perf_m(2,:)./sum(perf_m);
    
    figure('Position', [100 100 500 600], 'Color',[1 1 1]);
    hold on%subplot(2,1,1);
    for md = [1:numel(Poss_delays)-1]
        plot(Poss_dur, y_perf(:,md),'o','MarkerEdgeColor','none','MarkerFaceColor', [1 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))])
    end
    plot(Poss_dur, y_perf(:,end),'o','MarkerEdgeColor','none','MarkerFaceColor', [0 0 1])
    
    title('Behavioral performance on the iconic memory task');
    xlabel('stimulus duration (ms)')
    
    
    
%     clear Xreg flin Yreg
%     gestimates = 0.2;
%     for md = [1:numel(Poss_delays)]
%         Xreg{md} = Poss_dur;
%         str = ['@(param,xval)  param(' num2str(md+1) ')+param(1)*xval'];
%         flin_m = str2func(str);
%         flin{md} = flin_m;
%         Yreg{md} = AUC(md,:);
%         gestimates = [gestimates 20];
%     end
%     
    
    clear Xreg flin Yreg
    gestimates = [0.8 0.005];
    for md =  [1:numel(Poss_delays)]
        Xreg{md} = Poss_dur;
        str = ['@(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(' num2str(md+2) ')-xval)*param(2)))'];
        fsigm_m = str2func(str);
        fsigm{md} = fsigm_m;
        Yreg{md} = y_perf(:,md);
        gestimates = [gestimates 150];
    end
    
    [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(Xreg, Yreg, fsigm, gestimates);
    x1 = 0:10:Poss_dur(end);
    for md = [1:numel(Poss_delays)]-1
        y1 = 0.5+(beta(1)-0.5)./(1+10.^((beta(md+2)-x1)*beta(2)));
        plot(x1, y1, 'Color',[1 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))])
    end
    y1 = 0.5+(beta(1)-0.5)./(1+10.^((beta(end)-x1)*beta(2)));
    plot(x1, y1,'Color', [0 0 1])
    beta_behav = beta(3:end);
    
%     RSS2 = sum(r.^2);
    
    
%     % plot with shared fit parameters
%     subplot(2,1,2);plot(Poss_dur, y_perf)
%     if plot_mask
%         hold on;  plot(Poss_dur, y_perf_m,'r')
%     end
%     
%     fsigm = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(2)-xval)*param(3)));
%     fsigm_m = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(4)-xval)*param(3)));
%     [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_dur, Poss_dur}, {perf(2,:)./sum(perf),perf_m(2,:)./sum(perf_m)}, {fsigm, fsigm_m}, [.8  150 0.005 0.005]);
%     RSS1 = sum(r.^2);
%     F = ((RSS1-RSS2)/(6-4))/(RSS2/(numel(r)-6));
%     p = fcdf(F, 6-4, (numel(r)-6));
%     x1 = 0:10:Poss_dur(end);
%     y1 = 0.5+(beta(1)-0.5)./(1+10.^((beta(2)-x1)*beta(3)));
%     y1_m = 0.5+(beta(1)-0.5)./(1+10.^((beta(4)-x1)*beta(3)));
%     plot(x1, y1, 'Color','b')
%     plot(x1, y1_m, 'Color','r')
%     xlabel('stimulus duration (ms)')
    
    
    [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked', 'Location', 'SouthEast');
    LegendPos = get(legend_h,'Position');
    
%     set(legend_h, 'Location', 'NorthWest')
    set(legend_h, 'visible','off')
    EstimWorthBehav = (beta(4)-beta(2));
    dim = [LegendPos(1)-LegendPos(3) LegendPos(2) .3 LegendPos(4)/2];
    str = ['Estimated worth: ' num2str(EstimWorthBehav,3) 'ms'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
%     dim_F = dim; dim_F(2) = dim_F(2)+2*LegendPos(4); dim_F(4) = LegendPos(4);
%     if p<0.05
%         str = {'The more complex model explains the data', ['significantly better  (p = ' num2str(p,3) ')']};
%     else
%         str = {'The more complex model does', 'not explain the data',  ['significantly better (p = ' num2str(p,3) ')']};
%     end
%     annotation('textbox',dim_F,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
    
    filename = ['performances_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    %%
    f1 = figure('Position', [100 50 1400 900]);
    if plot_mask; f2 = figure('Position', [100 50 1400 900]);end
    for md = [1:numel(Poss_delays)]
    for sd = 1:numel(Poss_dur)
        t_connected_avg{sd,md} = 0*background{48,1,sd,md};n_t_connected_avg(sd)=0;
        t_unconnected_avg{sd,md} = 0*background{48,1,sd,md};n_t_unconnected_avg(sd)=0;
        for h = 1:2
            background_avg{h,sd,md} = 0*background{48,1,sd,md};n_background_avg(h,sd)=0;
            ttree_t_avg{h,sd,md} = 0*background{48,1,sd,md};n_ttree_t_avg(h,sd)=0;
            ttree_d_avg{h,sd,md} = 0*background{48,1,sd,md};n_ttree_d_avg(h,sd)=0;
            dtree_t_avg{h,sd,md} = 0*background{48,1,sd,md};n_dtree_t_avg(h,sd)=0;
            dtree_d_avg{h,sd,md} = 0*background{48,1,sd,md};n_dtree_d_avg(h,sd)=0;
        end
%         if plot_mask
%             t_connected_m_avg{sd,md} = 0*background{48,1,sd};n_t_connected_m_avg(sd)=0;
%             t_unconnected_m_avg{sd,md} = 0*background{48,1,sd};n_t_unconnected_m_avg(sd)=0;
%             for h = 2%1:
%                 background_m_avg{h,sd} = 0*background{48,1,sd};n_background_m_avg(h,sd)=0;
%                 ttree_t_m_avg{h,sd} = 0*background{48,1,sd};n_ttree_t_m_avg(h,sd)=0;
%                 ttree_d_m_avg{h,sd} = 0*background{48,1,sd};n_ttree_d_m_avg(h,sd)=0;
%                 dtree_t_m_avg{h,sd} = 0*background{48,1,sd};n_dtree_t_m_avg(h,sd)=0;
%                 dtree_d_m_avg{h,sd} = 0*background{48,1,sd};n_dtree_d_m_avg(h,sd)=0;
%             end
%         end
        for n = good_chann
            
            try t_connected_avg{sd,md} = nansum([t_connected_avg{sd,md} ,t_connected{n,sd,md}],2);end; if ~isnan(sum(t_connected{n,sd,md}));n_t_connected_avg(sd) = n_t_connected_avg(sd)+1;end
            try t_unconnected_avg{sd,md} =  nansum([t_unconnected_avg{sd,md},t_unconnected{n,sd,md}],2);end;if ~isnan(sum(t_unconnected{n,sd,md}));n_t_unconnected_avg(sd) = n_t_unconnected_avg(sd)+1;end
            for h = 1:2
                try background_avg{h,sd,md} =  nansum([background_avg{h,sd,md} ,background{n,h,sd,md}],2);end;if ~isnan(sum(background{n,h,sd,md}));n_background_avg(h,sd) = n_background_avg(h,sd)+1;end
                try ttree_t_avg{h,sd,md} =  nansum([ttree_t_avg{h,sd,md} ,ttree_t{n,h,sd,md}],2);end ;if ~isnan(sum(ttree_t{n,h,sd,md}));n_ttree_t_avg(h,sd) = n_ttree_t_avg(h,sd)+1;end
                try ttree_d_avg{h,sd,md} =  nansum([ttree_d_avg{h,sd,md} ,ttree_d{n,h,sd,md} ],2);end;if ~isnan(sum(ttree_d{n,h,sd,md}));n_ttree_d_avg(h,sd) = n_ttree_d_avg(h,sd)+1;end
                
                try dtree_t_avg{h,sd,md} =  nansum([dtree_t_avg{h,sd,md},dtree_t{n,h,sd,md}],2);end;if ~isnan(sum(dtree_t{n,h,sd,md}));n_dtree_t_avg(h,sd) = n_dtree_t_avg(h,sd)+1;end
                try dtree_d_avg{h,sd,md} =  nansum([dtree_d_avg{h,sd,md},dtree_d{n,h,sd,md}],2);end;if ~isnan(sum(dtree_d{n,h,sd,md}));n_dtree_d_avg(h,sd) = n_dtree_d_avg(h,sd)+1;end
            end
%             if plot_mask
%                 try t_connected_m_avg{sd,md} =  nansum([t_connected_m_avg{sd,md} ,t_connected_m{n,sd,md}],2);end;if ~isnan(sum(t_connected_m{n,sd,md}));n_t_connected_m_avg(sd) = n_t_connected_m_avg(sd)+1;end
%                 try t_unconnected_m_avg{sd,md} =  nansum([t_unconnected_m_avg{sd,md} ,t_unconnected_m{n,sd,md}],2);end;if ~isnan(sum(t_unconnected_m{n,sd,md}));n_t_unconnected_m_avg(sd) = n_t_unconnected_m_avg(sd)+1;end
%                 for h = 1:2
%                     try background_m_avg{h,sd} =  nansum([background_m_avg{h,sd} ,background_m{n,h,sd}],2);end;if ~isnan(sum(background_m{n,h,sd}));n_background_m_avg(h,sd) = n_background_m_avg(h,sd)+1;end
%                     try  ttree_t_m_avg{h,sd} =  nansum([ttree_t_m_avg{h,sd} ,ttree_t_m{n,h,sd}],2);end;if ~isnan(sum(ttree_t_m{n,h,sd}));n_ttree_t_m_avg(h,sd) = n_ttree_t_m_avg(h,sd)+1;end
%                     try ttree_d_m_avg{h,sd} =  nansum([ttree_d_m_avg{h,sd} ,ttree_d_m{n,h,sd} ],2);end;if ~isnan(sum(ttree_d_m{n,h,sd}));n_ttree_d_m_avg(h,sd) = n_ttree_d_m_avg(h,sd)+1;end
%                     
%                     try dtree_t_m_avg{h,sd} =  nansum([dtree_t_m_avg{h,sd},dtree_t_m{n,h,sd}],2);end;if ~isnan(sum(dtree_t_m{n,h,sd}));n_dtree_t_m_avg(h,sd) = n_dtree_t_m_avg(h,sd)+1;end
%                     try dtree_d_m_avg{h,sd} =  nansum([dtree_d_m_avg{h,sd},dtree_d_m{n,h,sd}],2);end;if ~isnan(sum(dtree_d_m{n,h,sd}));n_dtree_d_m_avg(h,sd) = n_dtree_d_m_avg(h,sd)+1;end
%                 end
%             end
        end
        
        t_connected_avg{sd,md} = t_connected_avg{sd,md}/n_t_connected_avg(sd);
        t_unconnected_avg{sd,md} = t_unconnected_avg{sd,md}/n_t_unconnected_avg(sd);
        for h = 1:2
            background_avg{h,sd,md} = background_avg{h,sd,md}/n_background_avg(h,sd);
            ttree_t_avg{h,sd,md} = ttree_t_avg{h,sd,md}/n_ttree_t_avg(h,sd);
            ttree_d_avg{h,sd,md} = ttree_d_avg{h,sd,md}/n_ttree_d_avg(h,sd);
            dtree_t_avg{h,sd,md} = dtree_t_avg{h,sd,md}/n_dtree_t_avg(h,sd);
            dtree_d_avg{h,sd,md} = dtree_d_avg{h,sd,md}/n_dtree_d_avg(h,sd);
        end
%         if plot_mask
%             t_connected_m_avg{sd,md} = t_connected_m_avg{sd,md}/n_t_connected_m_avg(sd);
%             t_unconnected_m_avg{sd,md} = t_unconnected_m_avg{sd,md}/n_t_unconnected_m_avg(sd);
%             for h = 1:2
%                 background_m_avg{h,sd,md} = background_m_avg{h,sd}/n_background_m_avg(h,sd);
%                 ttree_t_m_avg{h,sd} = ttree_t_m_avg{h,sd}/n_ttree_t_m_avg(h,sd);
%                 ttree_d_m_avg{h,sd} = ttree_d_m_avg{h,sd}/n_ttree_d_m_avg(h,sd);
%                 dtree_t_m_avg{h,sd} = dtree_t_m_avg{h,sd}/n_dtree_t_m_avg(h,sd);
%                 dtree_d_m_avg{h,sd} = dtree_d_m_avg{h,sd}/n_dtree_d_m_avg(h,sd);
%             end
%         end
        
        
        figure(f1);subplot(2,4,sd);hold on; title([num2str(Poss_dur(sd)), ' ms'])%
        %         plot(tb,smooth(background_avg{1,sd},20),'k');
        %         plot(tb,smooth(ttree_t_avg{1,sd},20),'b');
        %         plot(tb,smooth(ttree_d_avg{1,sd},20),'r');
        %         plot(tb,smooth(dtree_t_avg{1,sd},20),'c');
        %         plot(tb,smooth(dtree_d_avg{1,sd},20),'m');
        plot(tb,smooth(background_avg{2,sd,md},20),'k', 'LineWidth',2);
        plot(tb,smooth(ttree_t_avg{2,sd,md},20),'b', 'LineWidth',2);
        plot(tb,smooth(ttree_d_avg{2,sd,md},20),'r', 'LineWidth',2);
        plot(tb,smooth(dtree_t_avg{2,sd,md},20),'c', 'LineWidth',2);
        plot(tb,smooth(dtree_d_avg{2,sd,md},20),'m', 'LineWidth',2); xlim([-.2 0.9]);
        if sd == 1; legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'});end
        
%         if plot_mask
%             figure(f2);subplot(2,4,sd);hold on; title([num2str(Poss_dur(sd)), ' ms, masked at 20ms'])%
%             %             plot(tb,smooth(background_m_avg{1,sd},20),'k');
%             %             plot(tb,smooth(ttree_t_m_avg{1,sd},20),'b');
%             %             plot(tb,smooth(ttree_d_m_avg{1,sd},20),'r');
%             %             plot(tb,smooth(dtree_t_m_avg{1,sd},20),'c');
%             %             plot(tb,smooth(dtree_d_m_avg{1,sd},20),'m');
%             plot(tb,smooth(background_m_avg{2,sd},20),'k', 'LineWidth',2);
%             plot(tb,smooth(ttree_t_m_avg{2,sd},20),'b', 'LineWidth',2);
%             plot(tb,smooth(ttree_d_m_avg{2,sd},20),'r', 'LineWidth',2);
%             plot(tb,smooth(dtree_t_m_avg{2,sd},20),'c', 'LineWidth',2);
%             plot(tb,smooth(dtree_d_m_avg{2,sd},20),'m', 'LineWidth',2); xlim([-.2 0.9]);
%             if sd == 1; legend({'backgnd', 'target', 'target tree unconnected', 'distractor', 'distractor unconnected'});end
%             
%             
%         end
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
    
    
%     figure(f2)
%     subplot(2,4,sd);hold on; title([num2str(Poss_dur(sd)), ' ms, chan ' num2str(n)])
%     
%     yyaxis left
%     plot(tb, smooth(t_connected{n,sd,md}, 20),'Color', [1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays)) 1], 'LineStyle', '-', 'Marker','none' );hold on
%     plot(tb, smooth(t_unconnected{n,sd,md}, 20),'Color', [1 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))] , 'LineStyle', '-', 'Marker','none');
%     yyaxis right
%     plot(tb, smooth(t_connected{n,sd,md}-t_unconnected{n,sd,md}, 20),'Color',[1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))], 'LineStyle', '-', 'Marker','none'); xlim([-.2 0.9]);hold on
%     diffcu =smooth(t_connected{n,sd,md}-t_unconnected{n,sd,md},20);
%     %
% 
%     p=patch([0.01 1e-3*Poss_dur(sd)+integ_dur 1e-3*Poss_dur(sd)+integ_dur 0.01],[-10 -10 20 20],'r');
%     set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
% 
%     yyaxis left
%     set(gca, 'YLim', [-1 1.2*max(smooth(t_connected{n,sd,md}, 20))])
%     yyaxis right
%     set(gca, 'YLim', [min(diffcu)-0.2*(max(diffcu)-min(diffcu)) max(diffcu)+3*(max(diffcu)-min(diffcu))])
%     
    
    
ylim1 = [-2 2];    
for sd = 1:min(numel(Poss_dur),7)
     subplot(2,4,sd);hold on; title(['Stim duration: ' num2str(Poss_dur(sd)), ' ms']); xlabel('time(s)'); ylabel('normalized activity')
    for md = [1:numel(Poss_delays)]
       
        yyaxis left
        plot(tb, smooth(t_connected_avg{sd,md}, 20),'Color', [1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays)) 1], 'LineStyle', '-', 'Marker','none' );hold on;  ylim(ylim1)
        plot(tb, smooth(t_unconnected_avg{sd,md}, 20),'Color', [1 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))], 'LineStyle', '-', 'Marker','none' );
        yyaxis right
        plot(tb, smooth(t_connected_avg{sd,md}-t_unconnected_avg{sd,md}, 20),'Color',[1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays)) 1-(md/numel(Poss_delays))], 'LineStyle', '-', 'Marker','none'); xlim([-.2 0.9]);hold on
        line([-0.5 .9],[ 0 0])
        ylim([-.6 3])
        xlim([-.2 .9])
        if md==1 & sd==1
            legend({'no mask connected', 'no mask unconnected', 'conn - unconn'})
        end
        %         if sd == 1; legend({'conn masked', 'unconn masked', 'diff masked', 'conn ', 'unconn', 'diff '}); end
    end
    p=patch([0.01 1e-3*max(Poss_dur)+integ_dur 1e-3*max(Poss_dur)+integ_dur 0.01],[-10 -10 20 20],'r');
    set(p,'FaceAlpha',0.2, 'EdgeColor', 'none' );
end
    
    subplot(2,4,8); hold on
    %     integ_dur = 0.2;
    smoothness = 30;
    AUC = zeros(numel(Poss_delays),numel(Poss_dur));
    for md = [1:numel(Poss_delays)]
        for sd = 1:numel(Poss_dur)
            t_relev = find(tb>(0.01) & tb<((1e-3*max(Poss_dur))+integ_dur));
            diff_m =abs(smooth(t_connected_avg{sd,md}, smoothness)-smooth(t_unconnected_avg{sd,md}, smoothness));
            diff_m = diff_m(t_relev);%diff_m(diff_m<0)=0;
            %                     i = find(diff_m<0, 1);diff_m(i:end)=0;
            AUC(md,sd) = sum(diff_m);
            %             AUC_m(sd) = sum(smooth(t_connected_m_avg{sd,md}(t_relev)-t_unconnected_m_avg{sd,md}(t_relev), smoothness));
        end
        plot(Poss_dur, AUC(md,:), 'Color', [1-md/numel(Poss_delays) 1-md/numel(Poss_delays) 1-md/numel(Poss_delays)])
    end
    xlabel('Stimulus duration(ms)')
    
    
    try
        clear Xreg flin Yreg
        gestimates = 0.2;
        for md = [1:numel(Poss_delays)]
            Xreg{md} = Poss_dur;
            str = ['@(param,xval)  param(' num2str(md+1) ')+param(1)*xval'];
            flin_m = str2func(str);
            flin{md} = flin_m;
            Yreg{md} = AUC(md,:);
            gestimates = [gestimates 20];
        end
        %                 @(param,xval) param(1)+param(2)*xval;
        
        [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(Xreg, Yreg, flin, gestimates);
        x2 = 0:10:Poss_dur(end);
        for md = [1:numel(Poss_delays)-1]
            y2 = beta(md+1)+beta(1)*x2;
            plot(x2, y2, 'Color',[1 1-(md/numel(Poss_delays))  1-(md/numel(Poss_delays))])
            EstimWorthNeuronal(md,n) = (beta(end)-beta(md+1))/beta(1);
        end
        y2 = beta(end)+beta(1)*x2;
        plot(x2, y2, 'Color',[0 0 1])
        EstimWorthNeuronal(md+1,n) = (beta(end)-beta(md+2))/beta(1);
        
        [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked', 'Location', 'NorthWest');
        LegendPos = get(legend_h,'Position');
        set(legend_h,'visible','off')
        beta_neur = (30-beta(2:end))./beta(1);
        
        
        dim = [LegendPos(1) LegendPos(2)-LegendPos(4) .3 LegendPos(4)/2];
        str = ['Estimated worth: ' num2str(EstimWorthNeuronal(1,n),3) 'ms'];
        annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
        
    end
    
    
    
    
    title({'Integral of the difference in activity', ['between stim onset and ' num2str(integ_dur*1000) 'ms after the end of the stimulus']})

    filename = ['Diff_ConnectedUnconn_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    
    try
    figure; 
    plot(Poss_delays(1:end-1), beta_behav(1:end-1)-beta_behav(end),'b');
    hold on 
    plot(Poss_delays(1:end-1), beta_neur(1:end-1)-beta_neur(end),'r');
    title({'Estimated behavioral and neuronal worth of iconic memory', 'in function of the mask delay'})
    legend('behavior', 'neuronal')
    ylabel('Estimated worth (ms)')
    xlabel('delay between stim offset and mask onset(ms)')
   filename = ['CompWorthNeuroBehav_' num2str(Sessions(1)) '_' num2str(Sessions(end))];
    set(gcf, 'PaperPositionMode', 'auto')
    print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
    end
   
%     %% Compute sum of difference in the "masked condition" for the whole stim duration for all stim durations
%     Stim_info = zeros(1,numel(Poss_dur));
%     for sd = 1:min(numel(Poss_dur),7)
%         t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
%         diff_m =smooth(t_connected_avg{sd,md}, smoothness)-smooth(t_unconnected_avg{sd,md}, smoothness);
%         diff_m = diff_m(t_relev);
%         Stim_info(sd) = sum(diff_m);
%     end
%     figure;plot(Poss_dur, Stim_info, 'Color', [.1 .1 .1]);
%     
%     Stim_info_m = zeros(1,numel(Poss_dur));
%     for sd = 1:min(numel(Poss_dur),7)
%         t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
%         diff_m =smooth(t_connected_m_avg{sd,md}, smoothness)-smooth(t_unconnected_m_avg{sd,md}, smoothness);
%         diff_m = diff_m(t_relev);
%         Stim_info_m(sd) = sum(diff_m);
%     end
%     
%     hold on;plot(Poss_dur, Stim_info_m, 'Color', [0 .9 0]);
%     
%     flin = @(param,xval) param(1)+param(2)*xval;
%     flin_m = @(param,xval) param(3)+param(2)*xval;
%     [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_dur, Poss_dur}, {Stim_info,Stim_info_m}, {flin, flin_m}, [20 20 0.2]);
%     x2 = 0:10:Poss_dur(end);
%     y2 = beta(1)+beta(2)*x2;
%     y2_m =  beta(3)+beta(2)*x2;
%     plot(x2, y2, 'Color','k')
%     plot(x2, y2_m, 'Color','g')
%     [legend_h,object_h,plot_h,text_strings] = legend('no mask', 'masked', 'Location', 'NorthWest');
%     LegendPos = get(legend_h,'Position');
%     EstimWorthNeuronal = (beta(1)-beta(3))/beta(2);
%     dim = [LegendPos(1) LegendPos(2)-LegendPos(4) .3 LegendPos(4)/2];
%     str = ['Estimated worth: ' num2str(EstimWorthNeuronal,3) 'ms'];
%     annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none')
%     
%     title('Information about the connections in function of stimulus duration')
%     xlabel('stimulus duration')
%     
%     filename = ['Info_stimDur' num2str(Sessions(1)) '_' num2str(Sessions(end))];
%     set(gcf, 'PaperPositionMode', 'auto')
%     print(gcf, '-r0', [SaveDir filename, '.png'], '-dpng');
%     %%
    % compute integral of diff between connected and unconn curves over the
    % interval 50-400
    %     t_relev = find(tb>.01 & tb<0.5);
    if 0
        AUC = zeros(1,numel(Poss_dur));
        for sd = 1:min(numel(Poss_dur),7)
            t_relev = find(tb>(Poss_dur(sd)/1e3) & tb<(Poss_dur(sd)/1e3 +.3));
            AUC(sd) = sum(smooth(t_connected_avg{sd,md}(t_relev)-t_unconnected_avg{sd,md}(t_relev),20));
        end
        figure; plot(Poss_dur, AUC)
        if plot_mask
            AUC_m = zeros(1,numel(Poss_dur));
            for sd = 1:min(numel(Poss_dur),7)
                t_relev = find(tb>(Poss_dur(sd)/1e3) & tb<(Poss_dur(sd)/1e3 +.3));
                AUC_m(sd) = sum(smooth(t_connected_m_avg{sd,md}(t_relev)-t_unconnected_m_avg{sd,md}(t_relev), 20));
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
            cum_int(t) = sum(t_connected_avg{sd,md}(t_relev)-t_unconnected_avg{sd,md}(t_relev));
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
