% analysis function for the worth of iconic memory estimate 
% This script can be used for the initial worth estimation experiment
% (s9-11). It does not produce figures but performs the bootstrap
% statistics for the measures of iconic worth corrected for onset response
% of the unconnected condition
% written by Catherine Wacongne (catherine.waco@gmail.com)

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

% options
nboostrap = 1000;
integ_dur = 0.3;
smoothness = 30;
%9-11 is the basic iconic worth design.

Sessions = 9:11;


% get the conditions
isdrawn4 = [];
days =[];
for s = Sessions%
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
            Poss_dur = unique(stim_dur);
            mask = 2;
            if s>=15
                mask = 1+(2*LOG.Masked(trials).*(LOG.Mask_delay(trials)<60)+ (1-LOG.Masked(trials))); %.*(LOG.Mask_delay(trials)>60);
            else
                mask = mask.*(LOG.Masked(trials)) +1;
            end
        else
            isdrawn4 = [isdrawn4 isdrawn4_b];
            connected = [connected connected_b];
            targ = [targ LOG.targ_num(trials)];
            Hit = [Hit LOG.Hit(trials)];
            stim_dur = [stim_dur LOG.Stim_dur(trials)];
            Poss_dur = [Poss_dur unique(stim_dur)];
            mask_d = 2;
            if s>=15
                mask = [mask 1+(2*LOG.Masked(trials).*(LOG.Mask_delay(trials)<60)+ (1-LOG.Masked(trials)))];%LOG.Masked(trials).*(LOG.Mask_delay(trials)>60))];
            else
                mask =  [mask  mask_d.*(LOG.Masked(trials))+1];
            end
        end
        
    end
end
Poss_dur = unique(Poss_dur);
%%
Boostrap = randi(numel(mask),1000,numel(mask));%8030


SNR = zeros(1,48);
AllTrialsAllChan = cell(1,48);
for n = 1:48
    % Extract the trials from all blocs and conditions for each channel
    
    AllTrials = [];
    
    for s = Sessions%
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
EstimWorthBehav                 = zeros(1,nboostrap);
EstimWorthNeuronalAllChan       = zeros(1,nboostrap);
EstimWorthNeuronalAllChan_Relu  = zeros(1,nboostrap);
EstimWorthNeuronalAllChan_abs   = zeros(1,nboostrap);
EstimWorthNeuronalAllChanCorrected      = zeros(1,nboostrap);
EstimWorthNeuronalAllChanCorrected_Relu = zeros(1,nboostrap);
EstimWorthNeuronalAllChanCorrected_abs  = zeros(1,nboostrap);
for bst = 1:nboostrap
    if mod(bst,10)==0
        disp(['Bootstrap iteration ' num2str(bst) '/' num2str(nboostrap)])
    end
    bstrTrials = Boostrap(bst,:);
    bstdrawn4 = isdrawn4(bstrTrials);
    bstHit = Hit(bstrTrials);
    bstStim_dur = stim_dur(bstrTrials);
    bstmask = mask(bstrTrials);
    bsttarg = targ(bstrTrials);
    bstconnected = connected(bstrTrials);
    %% average the trials per channel for the bootstrap set
    background      = cell(48,2,numel(Poss_dur));
    perf            = zeros(2,numel(Poss_dur));
    perf_m          = zeros(2,numel(Poss_dur));
    t_connected     = cell(2,numel(Poss_dur));
    t_connected_m   = cell(2,numel(Poss_dur));
    t_unconnected   = cell(2,numel(Poss_dur));
    t_unconnected_m = cell(2,numel(Poss_dur));
    
    EstimWorthNeuronal      = zeros(1,48);
    EstimWorthNeuronal_Relu = zeros(1,48);
    EstimWorthNeuronal_abs  = zeros(1,48);
    OffResp_chan            = zeros(1,48);
    OnRespDistr_chan        = zeros(1,48);
    for n=1:48
        AllTrials = AllTrialsAllChan{n}(:, bstrTrials);
        for h = 1:2
            for sd = 1:numel(Poss_dur)
                % no mask conditions
                ind = bstdrawn4==0 & bstHit==h & bstStim_dur == Poss_dur(sd) & bstmask ==1 ;%& stim_dur == Poss_dur(sd) & bstmask == m
                background{n,h,sd}  = nanmean(AllTrials(:,ind),2);
                
                ind = find( bstHit==h & bstStim_dur == Poss_dur(sd) & bstmask ==1);
                perf(h,sd) = numel(ind);
                ind = bstconnected==1 & bstStim_dur == Poss_dur(sd) & bstmask ==1;% & Hit==2);%
                t_connected{n,sd} =  nanmean(AllTrials(:,ind),2);
                ind = bstconnected==2 & bstStim_dur == Poss_dur(sd) & bstmask ==1;%& Hit==2);%
                t_unconnected{n,sd} = nanmean(AllTrials(:,ind),2);
                
                % mask condition
                ind = find( bstHit==h & bstStim_dur == Poss_dur(sd) & bstmask ==3);
                perf_m(h,sd) = numel(ind);
                ind = bstconnected==1 & bstStim_dur == Poss_dur(sd) & bstmask ==3;% & Hit==2);%
                t_connected_m{n,sd} =  nanmean(AllTrials(:,ind),2);
                ind = find(bstconnected==2 & bstStim_dur == Poss_dur(sd) & bstmask ==3);% & Hit==2);%
                t_unconnected_m{n,sd} = nanmean(AllTrials(:,ind),2);
            end
        end
        
        %% channel level analysis: single channel neuronal worth and offset responses
        AUC = zeros(1,numel(Poss_dur)); %integral of the diff
        AUC_Relu = zeros(1,numel(Poss_dur)); % integral of positive parts of the diff
        AUC_abs = AUC_Relu;% integral of absolute value of the diff
        OffResp = AUC;
        OnRespDistr = AUC;
        for sd = 1:numel(Poss_dur)
            t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
            diff = smooth(t_connected{n,sd}(t_relev)-t_unconnected{n,sd}(t_relev),20);
            AUC(sd) = sum(diff);
            AUC_Relu(sd) = sum(diff(diff>0));
            AUC_abs(sd) = sum(abs(diff));
            
            tbl = find(tb>((1e-3*Poss_dur(sd))+0.01) & tb<((1e-3*Poss_dur(sd)) +0.03));
            blact = smooth(t_connected{n,sd}(tbl)-t_unconnected{n,sd}(tbl),20);
            pre_off = mean(blact);
            tpeak = find(tb>((1e-3*Poss_dur(sd))+0.03) & tb<((1e-3*Poss_dur(sd)) +0.09));
            pact = smooth(t_connected{n,sd}(tpeak)-t_unconnected{n,sd}(tpeak),20);
            peakvalue = max(pact);
            OffResp(sd) = peakvalue - pre_off;
            OnRespDistr(sd) = max(smooth(t_unconnected{n,sd}(pf)));
            %             OffRatio(sd) = peakvalue/pre_off;
        end
        AUC_m = zeros(1,numel(Poss_dur));
        AUC_Relu_m = AUC_m;
        AUC_abs_m = AUC_m;
        for sd = 1:numel(Poss_dur)
            t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
            diff_m =smooth(t_connected_m{n,sd}(t_relev)-t_unconnected_m{n,sd}(t_relev), 20);
            AUC_m(sd) = sum(diff_m);
            AUC_Relu_m(sd) = sum(diff_m(diff_m>0));
            AUC_abs_m(sd) = sum(abs(diff_m));
        end
        
        try
            flin = @(param,xval) param(1)+param(2)*xval;
            flin_m = @(param,xval) param(3)+param(2)*xval;
            [beta] = nlinmultifit({Poss_dur, Poss_dur}, {AUC,AUC_m}, {flin, flin_m}, [20 20 0.2]);
            [beta_Relu] = nlinmultifit({Poss_dur, Poss_dur}, {AUC_Relu,AUC_Relu_m}, {flin, flin_m}, [20 20 0.2]);
            [beta_abs] = nlinmultifit({Poss_dur, Poss_dur}, {AUC_abs,AUC_abs_m}, {flin, flin_m}, [20 20 0.2]);
            x2 = 0:10:Poss_dur(end);
            y2 = beta(1)+beta(2)*x2;
            y2_m =  beta(3)+beta(2)*x2;
        end
        EstimWorthNeuronal(n) = (beta(1)-beta(3))/beta(2);
        EstimWorthNeuronal_Relu(n) = (beta_Relu(1)-beta_Relu(3))/beta_Relu(2);
        EstimWorthNeuronal_abs(n) = (beta_abs(1)-beta_abs(3))/beta_abs(2);
        OffResp_chan(n) = mean(OffResp(end-2:end));
        OnRespDistr_chan(n) = mean(OnRespDistr(end-2:end));
        %         OffRatio_chann(n) = mean(OffRatio(end-2:end));
        
        
        
    end
    
    %% all chan analysis
    if 0
        [R,P] = corrcoef(OffResp_chan(good_chann),EstimWorthNeuronal(good_chann));
        linear_fit = polyfit(OffResp_chan(good_chann),EstimWorthNeuronal(good_chann),1);
        EstimWorthNeuronalAllChanCorrected(bst) = linear_fit(2);
        linear_fit = polyfit(OffResp_chan(good_chann),EstimWorthNeuronal_Relu(good_chann),1);
        EstimWorthNeuronalAllChanCorrected_Relu(bst) = linear_fit(2);
        linear_fit = polyfit(OffResp_chan(good_chann),EstimWorthNeuronal_abs(good_chann),1);
        EstimWorthNeuronalAllChanCorrected_abs(bst) = linear_fit(2);
    else 
        [R,P] = corrcoef(OnRespDistr_chan(good_chann),EstimWorthNeuronal(good_chann));
        linear_fit = polyfit(OnRespDistr_chan(good_chann),EstimWorthNeuronal(good_chann),1);
        EstimWorthNeuronalAllChanCorrected(bst) = linear_fit(2);
        linear_fit = polyfit(OnRespDistr_chan(good_chann),EstimWorthNeuronal_Relu(good_chann),1);
        EstimWorthNeuronalAllChanCorrected_Relu(bst) = linear_fit(2);
        linear_fit = polyfit(OnRespDistr_chan(good_chann),EstimWorthNeuronal_abs(good_chann),1);
        EstimWorthNeuronalAllChanCorrected_abs(bst) = linear_fit(2);
    end
    
    
    y_perf = perf(2,:)./sum(perf);
    y_perf_m = perf_m(2,:)./sum(perf_m);
    
    % fit with shared fit parameters
    fsigm = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(2)-xval)*param(3)));
    fsigm_m = @(param,xval) 0.5+(param(1)-0.5)./(1+10.^((param(4)-xval)*param(3)));
    [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit({Poss_dur, Poss_dur}, {perf(2,:)./sum(perf),perf_m(2,:)./sum(perf_m)}, {fsigm, fsigm_m}, [.8  150 0.005 0.005]);
    
    EstimWorthBehav(bst) = (beta(4)-beta(2));
    
    
    t_connected_avg     = cell(1, numel(Poss_dur));
    t_unconnected_avg   = cell(1, numel(Poss_dur));
    t_connected_m_avg   = cell(1, numel(Poss_dur));
    t_unconnected_m_avg = cell(1, numel(Poss_dur));
    n_t_connected_avg   = zeros(1,numel(Poss_dur));
    n_t_unconnected_avg = zeros(1,numel(Poss_dur));
    n_t_connected_m_avg = zeros(1,numel(Poss_dur));
    n_t_unconnected_m_avg = zeros(1,numel(Poss_dur));
    % grand averages
    for sd = 1:numel(Poss_dur)
        t_connected_avg{sd} = 0*background{48,1,sd};
        t_unconnected_avg{sd} = 0*background{48,1,sd};
        
        t_connected_m_avg{sd} = 0*background{48,1,sd};
        t_unconnected_m_avg{sd} = 0*background{48,1,sd};
        
        for n = good_chann
            try t_connected_avg{sd} = nansum([t_connected_avg{sd} ,t_connected{n,sd}],2);end; if ~isnan(sum(t_connected{n,sd}));n_t_connected_avg(sd) = n_t_connected_avg(sd)+1;end
            try t_unconnected_avg{sd} =  nansum([t_unconnected_avg{sd},t_unconnected{n,sd}],2);end;if ~isnan(sum(t_unconnected{n,sd}));n_t_unconnected_avg(sd) = n_t_unconnected_avg(sd)+1;end
            try t_connected_m_avg{sd} =  nansum([t_connected_m_avg{sd} ,t_connected_m{n,sd}],2);end;if ~isnan(sum(t_connected_m{n,sd}));n_t_connected_m_avg(sd) = n_t_connected_m_avg(sd)+1;end
            try t_unconnected_m_avg{sd} =  nansum([t_unconnected_m_avg{sd} ,t_unconnected_m{n,sd}],2);end;if ~isnan(sum(t_unconnected_m{n,sd}));n_t_unconnected_m_avg(sd) = n_t_unconnected_m_avg(sd)+1;end
        end
        
        t_connected_avg{sd} = t_connected_avg{sd}/n_t_connected_avg(sd);
        t_unconnected_avg{sd} = t_unconnected_avg{sd}/n_t_unconnected_avg(sd);
        t_connected_m_avg{sd} = t_connected_m_avg{sd}/n_t_connected_m_avg(sd);
        t_unconnected_m_avg{sd} = t_unconnected_m_avg{sd}/n_t_unconnected_m_avg(sd);
        
    end
    
    % Compute sum of difference in the "masked condition" for the whole stim duration for all stim durations
    Stim_info           = zeros(1,numel(Poss_dur));
    Stim_info_Relu      = zeros(1,numel(Poss_dur));
    Stim_info_abs       = zeros(1,numel(Poss_dur));
    Stim_info_Relu_m    = zeros(1,numel(Poss_dur));
    Stim_info_abs_m     = zeros(1,numel(Poss_dur));
    for sd = 1:min(numel(Poss_dur),7)
        t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
        diff_m =smooth(t_connected_avg{sd}, smoothness)-smooth(t_unconnected_avg{sd}, smoothness);
        diff_m = diff_m(t_relev);
        Stim_info(sd) = sum(diff_m);
        Stim_info_Relu(sd) = sum(diff_m(diff_m>0));
        Stim_info_abs(sd) = sum(abs(diff_m));
    end
    
    
    Stim_info_m = zeros(1,numel(Poss_dur));
    for sd = 1:min(numel(Poss_dur),7)
        t_relev = find(tb>(.01) & tb<((1e-3*Poss_dur(sd)) +integ_dur));
        diff_m =smooth(t_connected_m_avg{sd}, smoothness)-smooth(t_unconnected_m_avg{sd}, smoothness);
        diff_m = diff_m(t_relev);
        Stim_info_m(sd) = sum(diff_m);
        Stim_info_Relu_m(sd) = sum(diff_m(diff_m>0));
        Stim_info_abs_m(sd) = sum(abs(diff_m));
    end
    flin = @(param,xval) param(1)+param(2)*xval;
    flin_m = @(param,xval) param(3)+param(2)*xval;
    beta = nlinmultifit({Poss_dur, Poss_dur}, {Stim_info,Stim_info_m}, {flin, flin_m}, [20 20 0.2]);
    [beta_Relu] = nlinmultifit({Poss_dur, Poss_dur}, {Stim_info_Relu,Stim_info_Relu_m}, {flin, flin_m}, [20 20 0.2]);
    [beta_abs] = nlinmultifit({Poss_dur, Poss_dur}, {Stim_info_abs,Stim_info_abs_m}, {flin, flin_m}, [20 20 0.2]);
    
    
    EstimWorthNeuronalAllChan(bst) = (beta(1)-beta(3))/beta(2);
    EstimWorthNeuronalAllChan_Relu(bst) = (beta_Relu(1)-beta_Relu(3))/beta_Relu(2);
    EstimWorthNeuronalAllChan_abs(bst) = (beta_abs(1)-beta_abs(3))/beta_abs(2);
    
end
%%
mean_estimates = [mean(EstimWorthBehav(1:bst)) mean(EstimWorthNeuronalAllChan(1:bst)) mean(EstimWorthNeuronalAllChan_Relu(1:bst)) mean(EstimWorthNeuronalAllChan_abs(1:bst))];
std_estimates  = [std(EstimWorthBehav(1:bst)) std(EstimWorthNeuronalAllChan(1:bst)) std(EstimWorthNeuronalAllChan_Relu(1:bst)) std(EstimWorthNeuronalAllChan_abs(1:bst))];
figure; hold on
bar(1,mean_estimates(1), 'Facecolor', [.8 0 0]);
b = bar(2:4,mean_estimates(2:4), 'Facecolor', [0 1 0]);
errorbar(1:4,mean_estimates,std_estimates,'k.', 'LineWidth', 2)
set(gca, 'XTick', 1:4, 'XTickLabel', {'Behavior', 'Neur', 'Neur ReLu', 'Neur Abs'})

% b.CData(1,:) = [.5 0 .5];

mean_correstimates = [mean(EstimWorthBehav(1:bst)) mean(EstimWorthNeuronalAllChanCorrected(1:bst)) mean(EstimWorthNeuronalAllChanCorrected_Relu(1:bst)) mean(EstimWorthNeuronalAllChanCorrected_abs(1:bst))];
std_correstimates  = [std(EstimWorthBehav(1:bst)) std(EstimWorthNeuronalAllChanCorrected(1:bst)) std(EstimWorthNeuronalAllChanCorrected_Relu(1:bst)) std(EstimWorthNeuronalAllChanCorrected_abs(1:bst))];
figure; hold on
 bar(1,mean_correstimates(1), 'Facecolor', [.8 0 0]);
 bar(2:4,mean_correstimates(2:4), 'Facecolor', [0 1 0]);
errorbar(1:4,mean_correstimates,std_correstimates,'k.', 'LineWidth', 2)
set(gca, 'XTick', 1:4, 'XTickLabel', {'Behavior', 'Neur', 'Neur ReLu', 'Neur Abs'})

filename = ['BootstratEstimates_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
save(filename, 'EstimWorth*')

 %% p value tests 
a =  sort(EstimWorthBehav(1:bst)- EstimWorthNeuronalAllChan(1:bst) );
p = find(a>=0,1);if isempty(p); p=0;else; p = 1- abs((500-p)/500);end
pval = p; 

a =  sort(EstimWorthBehav(1:bst)- EstimWorthNeuronalAllChan_Relu(1:bst) );
p = find(a>=0,1);if isempty(p); p=0;else; p = 1- abs((500-p)/500);end
pval_Relu = p; 

a =  sort(EstimWorthBehav(1:bst)- EstimWorthNeuronalAllChan_abs(1:bst) );
p = find(a>=0,1);if isempty(p); p=0;else; p = 1- abs((500-p)/500);end
pval_abs = p; 

a =  sort(EstimWorthBehav(1:bst)- EstimWorthNeuronalAllChanCorrected(1:bst) );
p = find(a>=0,1);if isempty(p); p=0;else; p = 1- abs((500-p)/500);end
pval_Corr = p; 

a =  sort(EstimWorthBehav(1:bst)- EstimWorthNeuronalAllChanCorrected_Relu(1:bst) );
p = find(a>=0,1);if isempty(p); p=0;else; p = 1- abs((500-p)/500);end
pval_CorrRelu = p; 

a =  sort(EstimWorthBehav(1:bst)- EstimWorthNeuronalAllChanCorrected_abs(1:bst));
p = find(a>=0,1);if isempty(p); p=0;else; p = 1- abs((500-p)/500);end
pval_CorrAbs = p; 

absdiff1 = abs(EstimWorthBehav(1:bst)- EstimWorthNeuronalAllChanCorrected_abs(1:bst));
absdiff2 = abs(EstimWorthBehav(1:bst)- EstimWorthNeuronalAllChanCorrected(1:bst));

filename = ['BootstratPestimates_s' num2str(Sessions(1)) '_' num2str(Sessions(end))];
save(filename, 'pval*')
% [~,pval] = ttest(EstimWorthBehav(1:bst),EstimWorthNeuronalAllChan(1:bst));
% [~,pval_Relu] = ttest(EstimWorthBehav(1:bst),EstimWorthNeuronalAllChan_Relu(1:bst));
% [~,pval_abs] = ttest(EstimWorthBehav(1:bst),EstimWorthNeuronalAllChan_abs(1:bst));
% [~,pval_Corr] = ttest(EstimWorthBehav(1:bst),EstimWorthNeuronalAllChanCorrected(1:bst));
% [~,pval_CorrRelu] = ttest(EstimWorthBehav(1:bst),EstimWorthNeuronalAllChanCorrected_Relu(1:bst));
% [~,pval_CorrAbs] = ttest(EstimWorthBehav(1:bst),EstimWorthNeuronalAllChanCorrected_abs(1:bst));