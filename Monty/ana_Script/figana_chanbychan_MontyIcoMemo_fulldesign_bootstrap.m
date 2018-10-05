% function for the boostrap statistics corresponding to the full design 
% written by Catherine Wacongne 2018 (catherine.wacongne@gmail.com)

close all
clearvars

%directory of extracted data and stim log data
[extractdir, sessions, rawdir, SaveDir] = infoDir_MontyIcoMemo;
info = Log_MontyIcoMemo;

% time base
load([extractdir, sessions{1} '\EVT_', info(1).Tankname,'_block-',num2str(info(1).goodblocs(1))])
SF = EVENT.strms(1).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);

%options
smoothness = 30;
integ_dur = 0.3;
nboostrap = 1000;


%33-47 sessions with full design
Sessions = 33:47;%

% get the conditions
isdrawn4 = [];
days =[];
for s = Sessions
    for bloc = 1:numel(info(s).goodblocs)
        load([rawdir,'ICOmemo_' info(s).Tankname '_B',num2str(info(s).goodblocs(bloc))])
        trials =find(LOG.target_presented>0);
        isdrawn4_b = zeros(1,numel(trials));
            connected_b = zeros(1,numel(trials));
            for i =1: numel(trials)
                a = find(LOG.drawn_skeleton{trials(i)}==5);
                if ~isempty(a)
                    isdrawn4_b(i) = a;
                    connected_b(i) = LOG.drawn_conn{trials(i)}(a);
                end
                
            end
            days = [days, s*ones(1,numel(trials))];
            if isempty(isdrawn4)
                isdrawn4 = isdrawn4_b;
                connected = connected_b;
                targ = LOG.targ_num(trials);
                Hit = LOG.Hit(trials);
                stim_dur = LOG.Stim_dur(trials);
                mask = 2;
                mask = mask.*(LOG.Masked(trials)) +1;
                mask_delay = LOG.Mask_delay(trials);
                mask_delay(mask==1) = 500;

            else
                isdrawn4 = [isdrawn4 isdrawn4_b];
                connected = [connected connected_b];
                targ = [targ LOG.targ_num(trials)];
                Hit = [Hit LOG.Hit(trials)];
                stim_dur = [stim_dur LOG.Stim_dur(trials)];
                mask_d = 2;
                mask =  [mask  mask_d.*(LOG.Masked(trials))+1];
                mask_delay = [mask_delay LOG.Mask_delay(trials)];
                mask_delay(mask==1) = 500;
                                
            end
    end
end
Poss_dur = unique(stim_dur);
Poss_delays = unique(mask_delay);


%%

Boostrap = randi(numel(mask),1000,numel(mask));%8030

SNR = zeros(1,48);
AllTrialsAllChan = cell(1,48);
for n = 1:48
    % Extract the trials from all blocs and conditions for each channel
    
    AllTrials = [];
    
    for s = Sessions%60;%28:29%10:23%10:21;%10:12%6:numel(sessions)
        for bloc = 1:numel(info(s).goodblocs)
            %Read in the extracted data form the mapped channel
            load([extractdir sessions{s},'\Xtract_',info(s).Tankname,'_Block-',num2str(info(s).goodblocs(bloc)),'_',num2str(n)]);
            e = Env{1};
            
            if isempty(AllTrials)
                AllTrials=e;  
            else
                AllTrials = [AllTrials e];                
            end
            
        end
    end
   
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
     AllTrialsAllChan{n} = AllTrials2;
end 
good_chann = find(SNR>1);
%% Bootstrap
EstimWorthBehav                 = zeros(numel(Poss_delays),nboostrap);
EstimWorthNeuronalAllChan       = zeros(numel(Poss_delays),nboostrap);
EstimWorthNeuronalAllChan_Relu  = zeros(numel(Poss_delays),nboostrap);
EstimWorthNeuronalAllChan_abs   = zeros(numel(Poss_delays),nboostrap);
EstimWorthNeuronalAllChanCorrected      = zeros(numel(Poss_delays),nboostrap);
EstimWorthNeuronalAllChanCorrected_Relu = zeros(numel(Poss_delays),nboostrap);
EstimWorthNeuronalAllChanCorrected_abs  = zeros(numel(Poss_delays),nboostrap);
disp('Starting bootstrap...')    
for bst = 1:nboostrap
    if mod(bst,10)==0
        disp(['Bootstrap iteration ' num2str(bst) '/' num2str(nboostrap)])
    end
    bstrTrials = Boostrap(bst,:);
    bstdrawn4 = isdrawn4(bstrTrials);
    bstHit = Hit(bstrTrials);
    bstStim_dur = stim_dur(bstrTrials);
    bstmask = mask(bstrTrials);
    bstmask_delay = mask_delay(bstrTrials);
    bsttarg = targ(bstrTrials);
    bstconnected = connected(bstrTrials);
    %% average the trials per channel for the bootstrap set
    background      = cell(48,2,numel(Poss_delays),numel(Poss_dur));
    perf            = zeros(2,numel(Poss_delays),numel(Poss_dur));
    t_connected     = cell(2,numel(Poss_delays),numel(Poss_dur)); 
    t_unconnected   = cell(2,numel(Poss_delays),numel(Poss_dur));
    
    EstimWorthNeuronal      = zeros(numel(Poss_delays),48);
    EstimWorthNeuronal_Relu = zeros(numel(Poss_delays),48);
    EstimWorthNeuronal_abs  = zeros(numel(Poss_delays),48);
    OffResp_chan            = zeros(1,48);
    
    isbad_chann = zeros(1,48);
    for n=1:48
        AllTrials = AllTrialsAllChan{n}(:, bstrTrials);
        % get the conditions
        for h = 1:2
            for md = 1:numel(Poss_delays)
                for sd = 1:numel(Poss_dur)
                    ind = find( bstHit==h & bstStim_dur == Poss_dur(sd) & bstmask_delay==Poss_delays(md));
                    perf(h,md,sd) = numel(ind);
                                              
                    ind = bstdrawn4==0 & bstHit==h & bstStim_dur == Poss_dur(sd) & bstmask_delay==Poss_delays(md)  ;%& stim_dur == Poss_dur(sd) & bstmask == m
                    background{n,h,md,sd}  = nanmean(AllTrials(:,ind),2);

                    ind = bstconnected==1 & bstStim_dur == Poss_dur(sd) & bstmask_delay==Poss_delays(md);% & Hit==2);%
                    t_connected{n,md,sd} =  nanmean(AllTrials(:,ind),2);
                    
                    ind = bstconnected==2 & bstStim_dur == Poss_dur(sd) & bstmask_delay==Poss_delays(md);%& Hit==2);%
                    t_unconnected{n,md,sd} = nanmean(AllTrials(:,ind),2);
                     
                end
            end
        end
        %% channel level analysis: single channel neuronal worth 
        % get the AUC for all durations and delays, then computes a worth
        % for each delay
        % computes the size of the offset response (1 per chann)
        AUC = zeros(numel(Poss_delays),numel(Poss_dur)); %integral of the diff
        AUC_Relu = zeros(numel(Poss_delays),numel(Poss_dur)); % integral of positive parts of the diff
        AUC_abs = AUC_Relu;% integral of absolute value of the diff
        OffResp = zeros(1,numel(Poss_dur));
        
        for sd = 1:numel(Poss_dur)
            t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
            for md = 1:numel(Poss_delays)
                diff = smooth(t_connected{n,md, sd}(t_relev)-t_unconnected{n,md,sd}(t_relev),20);
                AUC(md,sd) = sum(diff);
                AUC_Relu(md,sd) = sum(diff(diff>0));
                AUC_abs(md,sd) = sum(abs(diff));
            end
            
            tbl = find(tb>((1e-3*Poss_dur(sd))+0.01) & tb<((1e-3*Poss_dur(sd)) +0.03));
            blact = smooth(t_connected{n,end,sd}(tbl)-t_unconnected{n,end,sd}(tbl),20);
            pre_off = mean(blact);
            tpeak = find(tb>((1e-3*Poss_dur(sd))+0.03) & tb<((1e-3*Poss_dur(sd)) +0.09));
            pact = smooth(t_connected{n,end,sd}(tpeak)-t_unconnected{n,end,sd}(tpeak),20);
            peakvalue = max(pact);
            OffResp(sd) = peakvalue - pre_off;
            %             OffRatio(sd) = peakvalue/pre_off;
 
        end
        OffResp_chan(n) = mean(OffResp(end-2:end));
        
        % lin regression for channel worth estimate
        try
            clear Xreg flin Yreg
            gestimates = 0.2;
            for md = 1:numel(Poss_delays)
                Xreg{md} = Poss_dur;
                str = ['@(param,xval)  param(' num2str(md+1) ')+param(1)*xval'];
                flin_m = str2func(str);
                flin{md} = flin_m;
                Yreg{md} = AUC(md,:);
                Yreg_Relu{md} = AUC_Relu(md,:);
                Yreg_abs{md} = AUC_abs(md,:);
                gestimates = [gestimates 20];
            end
            %                 @(param,xval) param(1)+param(2)*xval;
            
            [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(Xreg, Yreg, flin, gestimates);
            [beta_Relu,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(Xreg, Yreg_Relu, flin, gestimates);
            [beta_abs,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(Xreg, Yreg_abs, flin, gestimates);
            
            for md = 1:numel(Poss_delays)-1
                EstimWorthNeuronal(md,n) = (beta(end)-beta(md+1))/beta(1);
                EstimWorthNeuronal_Relu(md,n) = (beta_Relu(end)-beta_Relu(md+1))/beta_Relu(1);
                EstimWorthNeuronal_abs(md,n) = (beta_abs(end)-beta_abs(md+1))/beta_abs(1);
            end
        catch
            isbad_chann(n) = 1;
        end
        
        
    end
    %% Across chann analysis
    good_chann = find(SNR>1 & isbad_chann==0);
    for md = 1:numel(Poss_delays)-1
        linear_fit = polyfit(OffResp_chan(good_chann),EstimWorthNeuronal(md,good_chann),1);
        EstimWorthNeuronalAllChanCorrected(md,bst) = linear_fit(2);
        linear_fit = polyfit(OffResp_chan(good_chann),EstimWorthNeuronal_Relu(md,good_chann),1);
        EstimWorthNeuronalAllChanCorrected_Relu(md,bst) = linear_fit(2);
        linear_fit = polyfit(OffResp_chan(good_chann),EstimWorthNeuronal_abs(md,good_chann),1);
        EstimWorthNeuronalAllChanCorrected_abs(md,bst) = linear_fit(2);  
    end
    y_perf = squeeze(perf(2,:,:)./sum(perf));
    clear Xreg flin Yreg
    gestimates = [0.8 0.005];
    for md =  1:numel(Poss_delays)
        Xreg{md} = Poss_dur;
        str = ['@(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(' num2str(md+2) ')-xval)*param(2)))'];
        fsigm_m = str2func(str);
        fsigm{md} = fsigm_m;
        Yreg{md} = y_perf(md,:);
        gestimates = [gestimates 150];
    end
    [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(Xreg, Yreg, fsigm, gestimates);
     beta_behav = beta(3:end);
     EstimWorthBehav(1:numel(Poss_delays)-1,bst) = beta_behav(1:end-1)-beta_behav(end);
   
    
end
filename = ['BootstratEstimates_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
save(filename, 'EstimWorth*')