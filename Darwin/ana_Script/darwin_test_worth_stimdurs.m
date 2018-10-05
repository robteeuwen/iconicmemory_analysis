

% test darwin worth paradigm 

% load event 
q = 'EVT_Darwin_20180731_Block-1';
load(['D:\Projects\IconicMemory\datasets\Darwin\extractdata\20180731\' q]);

load(['D:\Projects\IconicMemory\datasets\Darwin\logs\ICOmemo_Darwin_20180731_B1.mat']);

% stim bits = stimulus onset, microbits = mask onset, target bits =
% response cue

ntrls = length(LOG.correct);

% find masked trials 
masked = find(LOG.isMask);

% find stimulus onset for each valid trial
stimon = nan(1,ntrls);
micro = nan(1,ntrls);
targ = nan(1,ntrls);
for i=1:length(EVENT.strons.Word)-1
    cw = EVENT.strons.Word(1,i); % current word
    nw = EVENT.strons.Word(1,i+1); % next word
    qs = find(EVENT.strons.Stim>cw & EVENT.strons.Stim<nw);
    qm = find(EVENT.strons.Micr>cw & EVENT.strons.Micr<nw);
    qt = find(EVENT.strons.Targ>cw & EVENT.strons.Targ<nw);
    if qs
        stimon(i) = EVENT.strons.Stim(qs);
    end
    if qm
        micro(i) = EVENT.strons.Micr(qm);
    end
    if qt 
        targ(i) = EVENT.strons.Targ(qt);
    end
end

% masked trials
stimdurs = micro(masked)-stimon(masked);
qstimdurs = LOG.Stim_dur(masked);
stimdurs = stimdurs*1000;

stimdurdif = qstimdurs - stimdurs;

stimdurs(isnan(stimdurs)) = [];


figure
histogram(stimdurs,30)

