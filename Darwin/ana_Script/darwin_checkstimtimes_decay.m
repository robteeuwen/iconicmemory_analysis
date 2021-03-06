% check stimulus times for darwin

close all
clearvars
dbstop if error

%% settings
rawdir = 'D:\Projects\IconicMemory\datasets\Darwin\RawData\';
datadir = 'D:\Projects\IconicMemory\datasets\Darwin\extractdata\';
logdir = 'D:\Projects\IconicMemory\datasets\Darwin\logs\';

% options
debug = 0;
refreshrate = 60; % refreshrate of the stimulus screen in Hz
date = '20181004';
block = 1;
baseline_period = 200; % firs ... ms, for photodiode drift correction
correctbaseline = 1;

% luminance levels of the diodesquare and their corresponding photodiode
% values
steplevels = [1 0.5 0.7 0]; % don't change this order (unless you know what you're doing)
lumlevels = [...
0,      0.001,      -0.01;...
0.5,    -0.05,  -0.075;...
0.7,    -0.085,  -0.15;...
1,      -0.37,  -0.5;];

lumlevels = [...
0,      0.001,      -0.1;...
0.5,    -0.15,  -0.3;...
0.7,    -0.31,  -0.6;...
1,      -0.7,  -2;];

% the TDT project records 4 analog inputs simultaneously, the diode should
% be plugged into channel 1
chan = 1;

datafile = ['XtractDio_Darwin_' date '_Block-' num2str(block)];
eventfile = ['EVT_Darwin_' date '_Block-' num2str(block)];

log = ['ICOmemo_Darwin_' date '_B' num2str(block)];

% import functions
addpath('C:\Users\teeuwen\Documents\MATLAB\standardfunctions');

%% load stuff
% load log
load([logdir log]);

% load extracted data
load([datadir date '\' datafile]);

% load event
load([datadir date '\' eventfile]);

% get x-axis for diode 
dioidx = find(strcmp('ANA1',{EVENT.strms.name}));
sampf = EVENT.strms(dioidx).sampf;
x = linspace(EVENT.Start,EVENT.Start+EVENT.Triallngth,EVENT.Triallngth*sampf);

d = Dio; %Dio{chan};

trials = find(LOG.target_presented>0);
ntrls = size(d,2);

%% get bits 
targbit = EVENT.strons.Targ;
stimbit = nan(1,length(targbit));

for t = 1:ntrls
    q = find(EVENT.strons.Stim<targbit(t),1,'last');
    stimbit(t) = EVENT.strons.Stim(q);
end

%% extract info from LOG
asked_stimdur = LOG.STIMDUR(trials);
asked_masktime = LOG.MASK_TIME(trials);


%% check refresh rate
powerspec = [];
for trl = 1:ntrls
    [powerspec(trl,:),f] = ez_powerspec(d(:,trl),sampf);
end

figure
range = find(f<100);
p = mean(powerspec,1);
plot(f(range),p(range));
title('Photodiode average powerspectrum over all trials');
xlabel('frequency(Hz)');

figure
plot(d)

%% convert diode timecourse to step function


% first get histogram of all luminances over all trials
base = (baseline_period/1000)*sampf;
allfp2 = cell(1,ntrls);
allfp = [];
for trl = 1:ntrls
    disp(['extracting peaks form trial ' num2str(trl)]);
    [fp, fpx] = findpeaks(d(:,trl)*-1,'MinPeakDistance',350);
    fp(1) = [];
    fpx(1) = [];
    
    fp = fp.*-1;
    
    baseline(trl) = mean(fp(find(fpx<base)));
    maxmin(trl) = max(fp)-min(fp);
    allfp2{trl} = fp;
    allfp = [allfp; fp];
end

baseline = smooth(baseline,round(ntrls/10));
baseline = baseline-mean(baseline);

if correctbaseline
    % correct for baseline
    allfp3 = [];
    for trl = 1:ntrls
        allfp3 = [allfp3; allfp2{trl}-baseline(trl)];
    end
    allfp = allfp3;
end

divs = 1;
lums = [];
for div = 1:divs
    qstart = (div-1)*(length(allfp)/divs)+1;
    qend = qstart+(length(allfp)/divs)-1;
    [n, bins] = histcounts(allfp(qstart:qend),200);
    cutoff = [];
    nstart = 1;
    while nstart < length(n)
        if sum(n(nstart:end)==0) > 0
            cutoff(end+1,1) = nstart+find(n(nstart:end)==0,1,'first')-1;
            nstart = cutoff(end,1);

            cutoff(end,2) = nstart+find(n(nstart:end)~=0,1,'first')-1;
            nstart = cutoff(end,2);
        else
            break
        end
    end
    edges = cutoff(:,1)+(cutoff(:,2)-cutoff(:,1))./2;
    lums(div,:) = sort(bins(round(edges)));
end
if size(lums,1)>1
    lums = mean(lums); 
end

figure
histogram(allfp,200); hold all;
set(gca,'YScale','log');
title('Shall we use these luminance levels as cutoffs?');
y_lim = get(gca,'YLim');
for l = lums
    plot([l l], y_lim);
end
pause

frametime = 1/refreshrate;
sampperframe = frametime*sampf;
badtrials = [];
dio_stimon = []; dio_stimoff = []; dio_maskon = []; dio_maskoff = []; dio_stimdur = []; dio_maskdur = []; maskdiff = []; stimdiff = [];
for trl = 1:ntrls

    disp(['analysing trial ' num2str(trl)]);
    c = d(:,trl);
    
    % find all peaks
    [fp, fpx] = findpeaks(c*-1,'MinPeakDistance',350);  
    fp = fp.*-1;
    fp = fp-baseline(trl);
    
    % remove the first peak
    fp(1) = [];
    fpx(1) = [];
    
    steps = nan(1,length(fpx));
    for ll = 0:length(lums)
        if ll == 0
            steps(find(fp<lums(1))) = steplevels(1);
        elseif ll == length(lums)
            steps(find(fp>lums(ll))) = steplevels(end); 
        else
            steps(find(fp>lums(ll) & fp<lums(ll+1))) = steplevels(ll+1);
        end
    end
    
    % make stepfunction using lumlevels - old approach
%     steps = nan(1,length(fpx));
%     for ll = 1:size(lumlevels,1)
%         steps(find(fp<lumlevels(ll,2) & fp>lumlevels(ll,3))) = lumlevels(ll,1);
%     end
    
    % make boxcar 
    boxc = [x(fpx(1)), steps(1)];
    for q = 2:length(steps)
        if steps(q) ~= boxc(end,2)
            boxc(end+1,:) = [x(fpx(q)) steps(q-1)];
            boxc(end+1,:) = [x(fpx(q)) steps(q)];
        end
    end
    boxc(end+1,:) = [x(fpx(end)) steps(end)];
    
    allbox{trl} = boxc;
    
    if sum(isnan(steps)) || debug
        figure; set(gcf,'Position',[234 107 1615 978]);
        plot(x,c); hold all; plot(x(fpx),fp,'ro');
        for ll = 1:length(lums)
            plot([x(fpx(1)) x(fpx(end))],[lums(ll) lums(ll)]);
        end
        disp(num2str(trl));
        warning('one or more of the photodiode values is outside the luminance levels boundaries you defined');
        pause
        close(gcf)
    end
    
    % find events 
    dio_stimon(trl) = x(fpx(find(steps==1,1,'first')));
    dio_stimoff(trl) = x(fpx(find(steps==1,1,'last')+1));
    dio_maskon(trl) = x(fpx(find(steps==0.5,1,'first')));
    dio_maskoff(trl) = x(fpx(find(steps==0.5,1,'last')+1));
    
    
    % CHECK STIMDUR and MASK ISI
    dio_stimdur(trl) = (dio_stimoff(trl)-dio_stimon(trl))*1000;
    stimdiff(trl) = abs(asked_stimdur(trl)-dio_stimdur(trl));
    
    dio_maskisi(trl) = (dio_maskon(trl)-dio_stimoff(trl))*1000;
    maskdiff(trl) = abs(asked_masktime(trl)-dio_maskisi(trl));
    if maskdiff(trl) > 2 || stimdiff(trl) > 2
        badtrials(end+1) = trl;
    end
end

% distribution for actual stimulus onset relative to stimbit 
figure
histogram(dio_stimon*1000);
xlabel('time difference (ms)');

% compare stimduration 
figure
subplot(1,2,1)
histogram(stimdiff);
subplot(1,2,2)
histogram(maskdiff);

% compute which trials are ok
goodtrials = ones(1,ntrls);
goodtrials(badtrials) = 0;

% plot boxcars 
figure
subplot(1,2,1); hold all; subplot(1,2,2); hold all;
cnd = unique(LOG.MASK_TIME);
cols = robcolors('moonrisekingdom',length(cnd));
for trl = 1:ntrls
    if goodtrials(trl)
        subplot(1,2,1)
    else
        subplot(1,2,2)
    end
    c = asked_masktime(trl);
    plot(allbox{trl}(:,1),allbox{trl}(:,2),'Color',cols(find(cnd==c),:));
end

% store in LOG
LOG.diode.goodtrials = goodtrials;
LOG.diode.dio_stimdur = dio_stimdur;
LOG.diode.dio_maskisi = dio_maskisi;
LOG.diode.stimdiff = stimdiff;
LOG.diode.maskdiff = maskdiff;
LOG.diode.stim_onset = dio_stimon;
save([logdir log],'LOG');

disp([num2str(sum(goodtrials)) ' trials out of ' num2str(length(goodtrials)) ' were ok.' ]);







