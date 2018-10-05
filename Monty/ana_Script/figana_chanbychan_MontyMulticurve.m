session = {'20171201'};

load(['D:\PostDoc\NeuronalData\Iconic\Monty\logs\Multicurve_Monty' session{1} '_B1.mat'])

load(['D:\PostDoc\NeuronalData\Iconic\Monty\extractdata\' session{1} '\EVT_Monty' session{1} '_Block-1.mat'])
SF = EVENT.strms(1).sampf;
TL = EVENT.Triallngth;
Start = EVENT.Start;
tb = (0:(TL*SF))./SF;
tb = tb+Start;
tb = tb(1:end-1);




trials = find(LOG.target_presented);
TargNum = LOG.targ_num(trials);
sideX = LOG.RFx(trials,:)>0;
sideY = LOG.RFy(trials,:)<0;
isinRF = (sideX+sideY)==2;
[targinRF,b] = find(isinRF');
RFonTargCurve = (targinRF'==TargNum);

good_chan = [27:30 33:44];
avg_conn = zeros(1,1373);
avg_unconn = zeros(1,1373);
for array = 1:2
    figure; 
    for n = (array-1)*24+1:array*24
        load(['D:\PostDoc\NeuronalData\Iconic\Monty\extractdata\' session{1} '\Xtract_Monty' session{1} '_Block-1_' num2str(n) '.mat'])
        e = Env{1};
        e = e';
        connected = nanmean(e(RFonTargCurve,:));
        unconnected = nanmean(e(RFonTargCurve==0,:));
        subplot(4,6,n-(array-1)*24);plot(tb,smooth(connected,20));hold on; plot(tb,smooth(unconnected,20),'r')
        xlim([-.2 1])
        if ~isempty(find(good_chan==n))
            avg_conn = avg_conn+connected;
            avg_unconn = avg_unconn+unconnected;
        end
    end
    
end

figure; plot(tb,smooth(avg_conn,20));hold on; plot(tb,smooth(avg_unconn,20),'r')
figure; plot(tb,smooth(avg_conn-avg_unconn,20), 'g');
