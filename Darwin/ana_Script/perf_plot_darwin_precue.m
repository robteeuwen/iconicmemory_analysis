% function perf_plot_darwin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that plots the performances of the monkey based on the log file
% writen by Catherine Wacongne - 03-05-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveDir = 'D:\PostDoc\NeuronalData\Iconic\Darwin\Figures\';
%% Log file location and meta info
rawdir = 'D:\PostDoc\NeuronalData\Iconic\Darwin\logs\'; % dir that contains the log files
info = Log_DarwinIcoMemo; % creates the logfile meta info (names and good blocs)
range_days = 94:95; % recording days that have to be taken into accont in the analysis

%% Reads the data and compute the performances

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
            target_num = LOG.targ_num(trials);
            stim_dur = LOG.Stim_dur(trials); % Stimulus duration
            Poss_dur = unique(stim_dur);
            CueTime = LOG.CueTime(trials);
            Poss_CueTimes = unique(CueTime);
            %             num_trees = cellfun(@numel,LOG.drawn_skeleton(trials)); % number of trees that were presented ( up to 5)
            
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
            target_num = [target_num LOG.targ_num(trials)];
            CueTime = [CueTime LOG.CueTime(trials)];
            Poss_CueTimes = [Poss_CueTimes unique(CueTime)];
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
            target_num = LOG.targ_num(trials);
            stim_dur = LOG.Stim_dur(trials); % Stimulus duration
            Poss_dur = unique(stim_dur);
            CueTime = LOG.CueTime(trials);
            Poss_CueTimes = unique(CueTime);
            %             num_trees = cellfun(@numel,LOG.drawn_skeleton(trials)); % number of trees that were presented ( up to 5)
            
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
            target_num = [target_num LOG.targ_num(trials)];
            CueTime = [ CueTime LOG.CueTime(trials)];
            Poss_CueTimes = [Poss_CueTimes unique(CueTime)];
            stim_dur = [stim_dur LOG.Stim_dur(trials)];
            Poss_dur = [Poss_dur unique(stim_dur)];
            %             num_trees = [num_trees  cellfun(@numel,LOG.drawn_skeleton(trials))];
        end
        Poss_dur = unique(Poss_dur);
        Poss_CueTimes = unique(Poss_CueTimes);
        
    end
    
    
end
Poss_dur = unique(Poss_dur);
Poss_CueTimes = unique(Poss_CueTimes);


 %% count the number of trials that belong to each condition 
 for h = 1:2 % hit
     for t_t = 1:5 % target tree        
         for sd = 1:numel(Poss_dur) % stim dur
             for ct = 1:numel(Poss_CueTimes)
                 ind = find( Hit==h & stim_dur == Poss_dur(sd) & target_tree==t_t & CueTime == Poss_CueTimes(ct) ); %target_tree==t_t &
                 perf(h,t_t,sd, ct) = numel(ind);
             end
             
         end         
     end
 end
 
 percent = squeeze(perf(2,:,:,:)./sum(perf));
 figure('Position', [100 50 1400 900]);
 for t = 1:5
     subplot(2,3,t); plot(Poss_CueTimes, squeeze(percent(t,1,:)),'r'); hold on; plot(Poss_CueTimes, squeeze(percent(t,2,:)), 'g');% plot(Poss_CueTimes, squeeze(percent(t,3,:)), 'b');
     ylim([.5 1])
     xlim([-50 120])
     xlabel('Cue Delay (ms)')
     ylabel('Performance')
     if t==1
         legend('short stim', 'long','Location', 'NorthEast')
     end
     title(['Perf when target belongs to Pair' num2str(t)]) 
 end



%% Cont trials per cond
perf = zeros(2,numel(Poss_dur), numel(Poss_CueTimes)); % correct(n/y) x target tree (1-5) x

for h = 1:2 % hit
    for sd = 1:numel(Poss_dur) % stim dur
        for ct = 1:numel(Poss_CueTimes)
            ind = find( Hit==h &  stim_dur == Poss_dur(sd) & CueTime == Poss_CueTimes(ct) ); %target_tree==t_t && num_trees==Poss_numTrees(n_t)
            perf(h,sd,ct) = numel(ind);
        end
    end
end

%% plot results
percent = squeeze(perf(2,:,:)./sum(perf));

perf_allDur = squeeze(sum(perf, 2));
percent_allDur = perf_allDur(2,:)./sum(perf_allDur);
figure; plot([-250 Poss_CueTimes(2:end)], percent_allDur); xlabel('Cue Time'); ylabel('performance');