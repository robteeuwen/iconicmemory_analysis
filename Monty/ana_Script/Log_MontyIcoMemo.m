function info = Log_MontyIcoMemo

% First log for Monty. Begining of training with stimuli tailored for both
% RF
i = 1;
info(i).Tankname = 'Monty20170216' ;
info(i).goodblocs = [5 6 7];

i = i+1;
info(i).Tankname = 'Monty20170428' ;
info(i).goodblocs = [1 2 3 4];

i = i+1;
info(i).Tankname = 'Monty20170503' ;
info(i).goodblocs = [1];

i = i+1;
info(i).Tankname = 'Monty20170504' ;
info(i).goodblocs = [1 4 5 6 7];

% good again
i = i+1;%5
info(i).Tankname = 'Monty20170511' ;
info(i).goodblocs = [1];

% switch to 3 pairs
i = i+1;
info(i).Tankname = 'Monty20170512' ;
info(i).goodblocs = [1];

i = i+1;
info(i).Tankname = 'Monty20170516' ;
info(i).goodblocs = [1];

i = i+1;
info(i).Tankname = 'Monty20170524' ;
info(i).goodblocs = [1 2];

% after a bit of training 
i = i+1;%9
info(i).Tankname = 'Monty20170531' ;
info(i).goodblocs = [1];

i = i+1;%10
info(i).Tankname = 'Monty20170601' ;
info(i).goodblocs = [2];

i = i+1;%11
info(i).Tankname = 'Monty20170602' ;
info(i).goodblocs = [1];


i = i+1; % attempt at varying mask SOA for stim =50ms bad perf 
info(i).Tankname = 'Monty20170606' ;
info(i).goodblocs = [1 2 ];

i = i+1; % mask at 100ms or 300ms after stim var stim length
info(i).Tankname = 'Monty20170607' ;
info(i).goodblocs = [1 ];


i = i+1; % mask at 60ms or 300ms after stim var stim length
info(i).Tankname = 'Monty20170608' ;
info(i).goodblocs = [1 ];

i = i+1; % mask at 60ms or 300ms after stim var stim length - timing shorter 
info(i).Tankname = 'Monty20170612' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20170613' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20170614' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20170615' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20170616' ;
info(i).goodblocs = [1 ];

i = i+1; % 80ms stim var mask delay
info(i).Tankname = 'Monty20170619' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20170620' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20170621' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20170622' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20170623' ;
info(i).goodblocs = [1 2];

i = i+1; % id
info(i).Tankname = 'Monty20170626' ;
info(i).goodblocs = [2];

i = i+1; % id
info(i).Tankname = 'Monty20170627' ;
info(i).goodblocs = [3];

i = i+1; % back to iconic worth test to see if mask still matters
info(i).Tankname = 'Monty20170704' ;
info(i).goodblocs = [1];

i = i+1; % iconic worth with mask at 0, 1frame, no msk WRONG
info(i).Tankname = 'Monty20170707' ;
info(i).goodblocs = [1];

i = i+1; % iconic worth with mask at 0, 1frame, no msk
info(i).Tankname = 'Monty20170710' ;
info(i).goodblocs = [4];

i = i+1; % iconic worth with mask at 0, 1frame, no msk
info(i).Tankname = 'Monty20170711' ;
info(i).goodblocs = [1];


i = i+1; % iconic worth with mask at various mask delay*stim dur (full design)
info(i).Tankname = 'Monty20170712' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170717' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170718' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170801' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170803' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170804' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170807' ;
info(i).goodblocs = [2];

i = i+1; % idem
info(i).Tankname = 'Monty20170808' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170809' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170810' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170811' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170814' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170816' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170817' ;
info(i).goodblocs = [2];

i = i+1; % idem
info(i).Tankname = 'Monty20170821' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170823' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170824' ;
info(i).goodblocs = [1];

i = i+1; % idem
info(i).Tankname = 'Monty20170829' ;
info(i).goodblocs = [1:4];

i = i+1; % idem
info(i).Tankname = 'Monty20170830' ;
info(i).goodblocs = [2];

i = i+1; % idem
info(i).Tankname = 'Monty20170831' ;
info(i).goodblocs = [1 2];

i = i+1; % idem
info(i).Tankname = 'Monty20170901' ;
info(i).goodblocs = [1 2];

i = i+1; % precue
info(i).Tankname = 'Monty20170904' ;
info(i).goodblocs = [1:5];

i = i+1; % precue 2 targ after 1 week of 1 targ training 
info(i).Tankname = 'Monty20170913' ;
info(i).goodblocs = [2];



i = i+1; % with 3 targ
info(i).Tankname = 'Monty20170922' ;
info(i).goodblocs = [1];


i = i+1; % with 2 targ
info(i).Tankname = 'Monty20170926' ;
info(i).goodblocs = [1];

i = i+1; % with 2 targ problem with the Stim B
info(i).Tankname = 'Monty20170927' ;
info(i).goodblocs = [1];

i = i+1; % with 2 targ
info(i).Tankname = 'Monty20170928' ;
info(i).goodblocs = [1:3];

i = i+1; % with 2 targ
info(i).Tankname = 'Monty20170929' ;
info(i).goodblocs = [1];

i = i+1; % with 2 targ
info(i).Tankname = 'Monty20171002' ;
info(i).goodblocs = [1];

i = i+1; % with 2 targ
info(i).Tankname = 'Monty20171003' ;
info(i).goodblocs = [1];

i = i+1; % with 2 targ
info(i).Tankname = 'Monty20171004' ;
info(i).goodblocs = [1 2];

i = i+1; % with 2 targ
info(i).Tankname = 'Monty20171005' ;
info(i).goodblocs = [1];

i = i+1; % with 2 targ
info(i).Tankname = 'Monty20171006' ;
info(i).goodblocs = [1];


i = i+1; % with 2 targ
info(i).Tankname = 'Monty20171016' ;
info(i).goodblocs = [3];

i = i+1; % with 2 targ corr timing
info(i).Tankname = 'Monty20171018' ;
info(i).goodblocs = [1];

i = i+1; % with 2 targ shorter stim dur and diff cue time 
info(i).Tankname = 'Monty20171019' ;
info(i).goodblocs = [1];

i = i+1; % with 3 targ shorter stim dur and diff cue time 
info(i).Tankname = 'Monty20171020' ;
info(i).goodblocs = [1:3];

i = i+1; % with 3 targ shorter stim dur and diff cue time 
info(i).Tankname = 'Monty20171023' ;
info(i).goodblocs = [1];


i = i+1; % id
info(i).Tankname = 'Monty20171024' ;
info(i).goodblocs = [1];

i = i+1; % changed cue from disappearance of non targets to brightening of candidate targets 
info(i).Tankname = 'Monty20171025' ;
info(i).goodblocs = [1 2];

i = i+1; % id
info(i).Tankname = 'Monty20171026' ;
info(i).goodblocs = [1];


i = i+1; % id
info(i).Tankname = 'Monty20171027' ;
info(i).goodblocs = [1];

i = i+1; % id
info(i).Tankname = 'Monty20171030' ;
info(i).goodblocs = [1];

i = i+1; % id
info(i).Tankname = 'Monty20171031' ;
info(i).goodblocs = [1];

i = i+1; % id
info(i).Tankname = 'Monty20171101' ;
info(i).goodblocs = [2 3];

i = i+1; % id
info(i).Tankname = 'Monty20171103' ;
info(i).goodblocs = [1:5];

i = i+1; % id with yellow distr
info(i).Tankname = 'Monty20171106' ;
info(i).goodblocs = [1];

i = i+1; % id with yellow distr
info(i).Tankname = 'Monty20171107' ;
info(i).goodblocs = [1];

i = i+1; % id with yellow distr
info(i).Tankname = 'Monty20171108' ;
info(i).goodblocs = [1 2];

i = i+1; % id with yellow distr
info(i).Tankname = 'Monty20171109' ;
info(i).goodblocs = [1 ];

i = i+1; % id with yellow distr const segm
info(i).Tankname = 'Monty20171110' ;
info(i).goodblocs = [1 ];

i = i+1; % id 
info(i).Tankname = 'Monty20171113' ;
info(i).goodblocs = [1 ];

i = i+1; % id 
info(i).Tankname = 'Monty20171114' ;
info(i).goodblocs = [1];

i = i+1; % id 
info(i).Tankname = 'Monty20171115' ;
info(i).goodblocs = [1:4];

i = i+1; % id 
info(i).Tankname = 'Monty20171116' ;
info(i).goodblocs = [1:3];

i = i+1; % id 
info(i).Tankname = 'Monty20171117' ;
info(i).goodblocs = [1 ];


i = i+1; % control normal CT bright cue 
info(i).Tankname = 'Monty20171121' ;
info(i).goodblocs = [1];


i = i+1; % control normal CT bright cue 
info(i).Tankname = 'Monty20171122' ;
info(i).goodblocs = [1 3 4 5];


i = i+1; % control normal CT bright cue 
info(i).Tankname = 'Monty20171123' ;
info(i).goodblocs = [1 2];


i = i+1; % control normal CT inverted contrats: black curves
info(i).Tankname = 'Monty20171124' ;
info(i).goodblocs = [1];

i = i+1; % control normal CT bright cue - biais towards pair1
info(i).Tankname = 'Monty20171127' ;
info(i).goodblocs = [1];

i = i+1; % control normal CT bright cue - biais towards pair5
info(i).Tankname = 'Monty20171128' ;
info(i).goodblocs = [1];

i = i+1; % control normal CT bright cue - biais towards pair5
info(i).Tankname = 'Monty20171129' ;
info(i).goodblocs = [1];

% i = i+1; % control runstim_multicurve (4 curves)
% info(i).Tankname = 'Monty20171130' ;
% info(i).goodblocs = [1];
% 
% i = i+1; % control runstim_multicurve (4 curves)
% info(i).Tankname = 'Monty201711201' ;
% info(i).goodblocs = [1];
% 
% i = i+1; % control runstim_multicurve (6 curves)
% info(i).Tankname = 'Monty20171204' ;
% info(i).goodblocs = [2];
% 
i = i+1; % new precue based on multicurve
info(i).Tankname = 'Monty20180112' ;
info(i).goodblocs = [2];

i = i+1; % new precue based on multicurve
info(i).Tankname = 'Monty20180115' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve
info(i).Tankname = 'Monty20180116' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve
info(i).Tankname = 'Monty20180117' ;
info(i).goodblocs = [1:5];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10
info(i).Tankname = 'Monty20180118' ;
info(i).goodblocs = [1:4];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10
info(i).Tankname = 'Monty20180119' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10
info(i).Tankname = 'Monty20180122' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10
info(i).Tankname = 'Monty20180123' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10
info(i).Tankname = 'Monty20180124' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10
info(i).Tankname = 'Monty20180125' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10
info(i).Tankname = 'Monty20180126' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10
info(i).Tankname = 'Monty20180129' ;
info(i).goodblocs = [1];

i = i+1; % new precue based on multicurve with cue delay -4 0 4 10 15 + longer disapearing segm
info(i).Tankname = 'Monty20180130' ;
info(i).goodblocs = [1];

i = i+1; % id
info(i).Tankname = 'Monty20180131' ;
info(i).goodblocs = [1 2];

i = i+1; % id
info(i).Tankname = 'Monty20180201' ;
info(i).goodblocs = [1 2];

i = i+1; % id
info(i).Tankname = 'Monty20180202' ;
info(i).goodblocs = [1 ];

i = i+1; % id
info(i).Tankname = 'Monty20180205' ;
info(i).goodblocs = [1];

i = i+1; % id
info(i).Tankname = 'Monty20180206' ;
info(i).goodblocs = [1];

i = i+1; % id
info(i).Tankname = 'Monty20180207' ;
info(i).goodblocs = [1 2];

i = i+1; % id
info(i).Tankname = 'Monty20180208' ;
info(i).goodblocs = [1];

i = i+1; % id
info(i).Tankname = 'Monty20180209' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180212' ;
info(i).goodblocs = [2];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180214' ;
info(i).goodblocs = [1:5];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180215' ;
info(i).goodblocs = [2];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180216' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180219' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180220' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180221' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180222' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs
info(i).Tankname = 'Monty20180223' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs variable stim dur, no mask
info(i).Tankname = 'Monty20180227' ;
info(i).goodblocs = [1];

% i = i+1; % disapp segm in RFs variable stim dur & mask
% info(i).Tankname = 'Monty20180228' ;
% info(i).goodblocs = [1];
% 
% i = i+1; % disapp segm in RFs variable stim dur & mask
% info(i).Tankname = 'Monty20180301' ;
% info(i).goodblocs = [1 2];
% 
% 
% i = i+1; % disapp segm in RFs variable stim dur & mask
% info(i).Tankname = 'Monty20180302' ;
% info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs variable stim dur & mask
info(i).Tankname = 'Monty20180305' ;
info(i).goodblocs = [1];


i = i+1; % disapp segm in RFs variable stim dur & mask
info(i).Tankname = 'Monty20180306' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs variable stim dur & mask
info(i).Tankname = 'Monty20180307' ;
info(i).goodblocs = [1:5];

i = i+1; % disapp segm in RFs variable stim dur & mask
info(i).Tankname = 'Monty20180308' ;
info(i).goodblocs = [1:6];

i = i+1; % disapp segm in RFs variable stim dur & mask
info(i).Tankname = 'Monty20180312' ;
info(i).goodblocs = [2];

i = i+1; % disapp segm in RFs variable stim dur & mask
info(i).Tankname = 'Monty20180313' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs variable stim dur & mask
info(i).Tankname = 'Monty20180314' ;
info(i).goodblocs = [1:5];

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180315' ;
info(i).goodblocs = [1:3];

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180316' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180319' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180320' ;
info(i).goodblocs = [1];

%% ROB STARTED EDITING THIS DOCUMENT FROM HERE ONWARDS ON THE 5TH OF APRIL 2018

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180321' ;
info(i).goodblocs = [1:4];

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180322' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180409' ;
info(i).goodblocs = [1:5];

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180410' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs fix stim dur & variable mask onset
info(i).Tankname = 'Monty20180411' ;
info(i).goodblocs = [1:5];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180412' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180413' ;
info(i).goodblocs = [1:5];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180416' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180417' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180418' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180419' ;
info(i).goodblocs = [1:3];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180420' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180423' ;
info(i).goodblocs = [1];

i = i+1; % disapp segm in RFs SHORTER stim dur & variable mask onset
info(i).Tankname = 'Monty20180424' ;
info(i).goodblocs = [1];