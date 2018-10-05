% function perf_plot_darwin 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that plots the performances of the monkey based on the log file 
% writen by Catherine Wacongne - 03-05-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveDir = 'C:\Users\wacongne\Documents\PostDoc\NeuronalData\Iconic\Darwin\Figures\';
%% Log file location and meta info
rawdir = 'C:\Users\wacongne\Documents\PostDoc\NeuronalData\Iconic\Darwin\logs\'; % dir that contains the log files
info = Log_DarwinIcoMemo; % creates the logfile meta info (names and good blocs)
range_days = 24:30%24:29;%19:23;%24:28;%23;%15:18;%5:9; % recording days that have to be taken into accont in the analysis 

%% Reads the data and compute the performances 
perf = zeros(2,3, 2,numel(range_days)); % correct(n/y) x target tree (1-5) x 
                        %mask type (no,cue,mask) x stim_duration x
                        %#targ(3-5)x days
AllCond = [];
 for s=1:numel(range_days) 
     %% Gets all the conditions for the good trials 
    try 
        load([rawdir,'ICOmemo_' info(range_days(s)).Tankname '_B0'])
        trials =find(LOG.target_presented>0);
        if isempty(AllCond) % for the 1st good bloc of the day
            AllCond=trials;           
            Hit = LOG.correct(trials)*2 + LOG.congruent_error(trials)+LOG.incongruent_error(trials) + LOG.incongruent_hit(trials);
            congruent_error = LOG.congruent_error(trials);
            incongruent_hit = LOG.incongruent_hit(trials);
            incongruent_error = LOG.incongruent_error(trials);
            dsk = LOG.drawn_skeleton(trials);
            target_tree = zeros(1,numel(trials)); % 1 to 5
%             numtrees = cellfun(@numel,LOG.drawn_skeleton(trials));
            for i=1:numel(trials)
                target_tree(i) = dsk{i}(LOG.chosen_Cue{trials(i)});
            end
            mask = LOG.ismask(trials)+1;
            mask = mask.*(LOG.CueDuration(trials)==100)-  (LOG.CueDuration(trials)==100).*(LOG.CueTime(trials)>10) + 1;
            mask_type = mask;%(LOG.ismask(trials)+1) .* (LOG.CueDuration(trials)==100) +1; % 1 = no mask; 2 = post cue; 3 = mask
            stim_dur = LOG.Stim_dur(trials); % Stimulus duration 
            Poss_dur = unique(stim_dur); 
            num_trees = cellfun(@numel,LOG.drawn_skeleton(trials)); % number of trees that were presented ( up to 5)
                        
        else % append info for the next blocs
            AllCond = [AllCond trials];
            Hit = [Hit (LOG.correct(trials)*2  + LOG.congruent_error(trials)+LOG.incongruent_error(trials) + LOG.incongruent_hit(trials))];

            congruent_error = [congruent_error LOG.congruent_error(trials)];
            incongruent_hit = [incongruent_hit LOG.incongruent_hit(trials)];
            incongruent_error = [incongruent_error LOG.incongruent_error(trials)];
            dsk = LOG.drawn_skeleton(trials);
            target_tree_tmp = zeros(1,numel(trials));
            for i=1:numel(trials)
                 target_tree_tmp(i) = dsk{i}(LOG.chosen_Cue{trials(i)});
            end

            target_tree = [target_tree target_tree_tmp];
            mask = LOG.ismask(trials)+1;
            mask = mask.*(LOG.CueDuration(trials)==100) -  (LOG.CueDuration(trials)==100).*(LOG.CueTime(trials)>10)+ 1;
            mask_type = [mask_type mask];%(LOG.ismask(trials)+1).*(LOG.CueDuration(trials)==100)+1];
            stim_dur = [stim_dur LOG.Stim_dur(trials)];
            Poss_dur = [Poss_dur unique(stim_dur)];
            num_trees = [num_trees  cellfun(@numel,LOG.drawn_skeleton(trials))];
        end
    end
    for bloc = 1:numel(info(range_days(s)).goodblocs)
        load([rawdir,'ICOmemo_' info(range_days(s)).Tankname '_B',num2str(info(range_days(s)).goodblocs(bloc))])
        trials =find(LOG.target_presented>0);
        if isempty(AllCond) % for the 1st good bloc of the day
            AllCond=trials;           
            Hit = LOG.correct(trials)*2 + LOG.congruent_error(trials)+LOG.incongruent_error(trials) + LOG.incongruent_hit(trials);
            congruent_error = LOG.congruent_error(trials);
            incongruent_hit = LOG.incongruent_hit(trials);
            incongruent_error = LOG.incongruent_error(trials);
            dsk = LOG.drawn_skeleton(trials);
            target_tree = zeros(1,numel(trials)); % 1 to 5
%             numtrees = cellfun(@numel,LOG.drawn_skeleton(trials));
            for i=1:numel(trials)
                target_tree(i) = dsk{i}(LOG.chosen_Cue{trials(i)});
            end
            mask = LOG.ismask(trials)+1;
            mask = mask.*(LOG.CueDuration(trials)==100)-  (LOG.CueDuration(trials)==100).*(LOG.CueTime(trials)>60) + 1;
            mask_type = mask;%(LOG.ismask(trials)+1) .* (LOG.CueDuration(trials)==100) +1; % 1 = no mask; 2 = post cue; 3 = mask
            stim_dur = LOG.Stim_dur(trials); % Stimulus duration 
            Poss_dur = unique(stim_dur); 
            num_trees = cellfun(@numel,LOG.drawn_skeleton(trials)); % number of trees that were presented ( up to 5)
                        
        else % append info for the next blocs
            AllCond = [AllCond trials];
            Hit = [Hit (LOG.correct(trials)*2  + LOG.congruent_error(trials)+LOG.incongruent_error(trials) + LOG.incongruent_hit(trials))];

            congruent_error = [congruent_error LOG.congruent_error(trials)];
            incongruent_hit = [incongruent_hit LOG.incongruent_hit(trials)];
            incongruent_error = [incongruent_error LOG.incongruent_error(trials)];
            dsk = LOG.drawn_skeleton(trials);
            target_tree_tmp = zeros(1,numel(trials));
            for i=1:numel(trials)
                 target_tree_tmp(i) = dsk{i}(LOG.chosen_Cue{trials(i)});
            end

            target_tree = [target_tree target_tree_tmp];
            mask = LOG.ismask(trials)+1;
            mask = mask.*(LOG.CueDuration(trials)==100) -  (LOG.CueDuration(trials)==100).*(LOG.CueTime(trials)>60)+ 1;
            mask_type = [mask_type mask];%(LOG.ismask(trials)+1).*(LOG.CueDuration(trials)==100)+1];
            stim_dur = [stim_dur LOG.Stim_dur(trials)];
            Poss_dur = [Poss_dur unique(stim_dur)];
            num_trees = [num_trees  cellfun(@numel,LOG.drawn_skeleton(trials))];
        end
        Poss_dur = unique(Poss_dur);
        Poss_numTrees = unique(num_trees);
        
    end
    %% count the number of trials that belong to each condition 
    for h = 1:2 % hit
%         for t_t = 1:5 % target tree
            for m_t = 1:3 % mask type
                for sd = 1:numel(Poss_dur) % stim dur
%                     for n_t = 1:numel(unique(num_trees)) % number of trees
                        ind = find( Hit==h &  mask_type == m_t & stim_dur == Poss_dur(sd)  ); %target_tree==t_t && num_trees==Poss_numTrees(n_t)
                        perf(h,m_t,sd,s) = numel(ind);
%                     end
                end
            end
%         end
    end
   
    
    
    
    
 end
 
 
 %% error analysis
    figure('Position', [100 50 1400 900]);
    for tt = 1:5
        for m = [1 3]
            subplot(5,2,(tt-1)*2 + 1 +(m==3));hold on
            corr = numel(find(Hit==2 & mask_type==m & target_tree==tt ));
            targt_err = numel(find(Hit==1 & mask_type==m & target_tree==tt & congruent_error==1));
            distt_conn =numel(find(Hit==1 & mask_type==m & target_tree==tt & incongruent_hit==1));
            distt_err = numel(find(Hit==1 & mask_type==m & target_tree==tt & incongruent_error==1));
            bar(1,corr,'g'); bar(2,targt_err,'r');bar(3,distt_conn,'c'); bar(4,distt_err,'m')
            
            set(gca,'XTick',1:4)
            set(gca,'XTickLabel',{'corr','c distr','uc targ','uc d'});
            ylabel([ {'cued pair '}; {num2str(tt)}])
            
            if tt==1
                if m==1
                    title('no mask')
                else
                    title('mask')
                end
            end
        end
    end
    
    set(gcf,'PaperPositionMode','auto')
    print(gcf,'-dpng',[SaveDir 'TargetChoice_s' num2str(range_days(1)) '_' num2str(range_days(end)) '.png'])

 % number of trials when not taking into account the target tree and number
 % of trees
% average accross days
perf = sum(perf,4); 
%% Makes the figures


 s_perf = squeeze(perf(2,:,:)./sum(perf));
 
    
    % plot of perf in function of sd for the 3 types of mask averaged
    % across all target tree and number of trees 
    figure;
%     sd_perf = squeeze(nanmean(s_perf,1));
    plot(Poss_dur, squeeze(s_perf(1,:))', Poss_dur, squeeze(s_perf(3,:)), 'r')%Poss_dur, squeeze(s_perf(2,:)), 'g'
    legend({'no mask', 'mask'})
    xlabel('stimulus duration')
    ylabel('performance')
    hold on
    [param]=sigm_fit([ Poss_dur],[squeeze(s_perf(1,:))], [.5, NaN, NaN, NaN], [.5 .8  50 0.005], 0);
    y=fsigm(param,0:10:Poss_dur(end));
    plot(0:10:Poss_dur(end),y, 'Color', 'b'); 
    
    Half_dur = param(3);
    param
    
    [param]=sigm_fit([ Poss_dur],[squeeze(s_perf(3,:))], [0.5, NaN, NaN, NaN], [.5 .8  50 0.005], 0);
    y=fsigm(param,0:10:Poss_dur(end));
    plot(0:10:Poss_dur(end),y, 'Color', 'r'); 
    
     Half_dur_m = param(3);
     param_m = param