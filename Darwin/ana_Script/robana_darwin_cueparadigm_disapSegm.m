
clear all

% get some info and load an EVENT structure
[extractdir, sessions, rawdir, SaveDir] = infoDir_DarwinIcoMemo;
info = Log_DarwinIcoMemo;
load([extractdir, sessions{end} '\EVT_', info(numel(sessions)).Tankname,'_block-',num2str(info(numel(sessions)).goodblocs(1))])

% get some more info from the EVENT structure
SF = EVENT.strms(1).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);
StimCondStr = {'Short', 'Long'};
Cond = cell(96,2,2,2);
ci = zeros(96,2);

% set sessions 
Sessions = 111:116;%
SNR = zeros(1,48);

% target connected in RF: 6
% target unconnected in RF: 5 (target pair distractor)
curveInRF = 6; % pair in RF, curve in RF
curveOutRF = 5; % pair in RF, curve out RF
ana_chan = 25:48;

% timing: 
% precue skeleton duration: 247 ms 
% connector duration: 294 ms
% cue onset: 
% - 0 ms after stimulus offset
% - 47 ms after stimulus offset
% - 117 ms after stimulus offset
% - 176 ms after stimulus offset

for n = ana_chan
    
    AllTrials = [];
    AllCond = [];
    days =[];
    for s = Sessions
        for bloc = 1:numel(info(s).goodblocs)
            
            % load the data 
            load([extractdir sessions{s},'\Xtract_',info(s).Tankname,'_Block-',num2str(info(s).goodblocs(bloc)),'_',num2str(n)]); % data
            load([extractdir, sessions{s} '\EVT_', info(s).Tankname,'_block-',num2str(info(s).goodblocs(bloc))]) % EVENT
            load([rawdir,'ICOmemo_' info(s).Tankname '_B',num2str(info(s).goodblocs(bloc))]) % LOG
            
            e = Env{1}; % envelope data in time x trials
            trials =find(LOG.target_presented>0);
            
            % correct number of trials
            if numel(trials)>size(e,2)
                trials = trials(1:size(e,2));
            end
            
            cuetimes = LOG.CueDelay(trials);
            allcuetimes = sort(unique(cuetimes));
            
            targ = LOG.targ_num(trials);
            islong = LOG.conn_codes(trials,1)'-1;
            
            for c = 1:length(allcuetimes)
                ct = allcuetimes(c); 
                ctrl = cuetimes==ct; 
                
                % cued tree in RF, connected branch in RF
                q = find(targ==CurveInRF && ctrl);
                
                % cued tree in Rf, unconnected branch in RF
                
                % uncued tree in Rf, connected branch in RF
                
                % uncued tree in RF, unconnected branch in RF
            end
        end
    end
end


