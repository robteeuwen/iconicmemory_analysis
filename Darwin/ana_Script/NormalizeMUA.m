function CleanData = NormalizeMUA(MUA, chan, days,hit ,bf, pf)

%Normalisation approach
%Lets try normliazing per day
MUA = MUA';
disp(['Normalizing the data of channel ',num2str(chan)])
n_days = unique(days);
for FL = n_days
    
    %get data from this day that isn;t an outlier
    %             f = find(days == FL & ~badtrials);
    
    %get data from this day
    f = find(days == FL);
    
    %NOW remove crazy stuff
    %Sample by sample outlier removal
    buf = MUA(f,:);
    MUAv = reshape(buf,size(buf,1)*size(buf,2),1);
    out = abs((MUAv-nanmean(MUAv))./nanstd(MUAv))>5;
    MUAv(out) = NaN;
    out = abs((MUAv-nanmean(MUAv))./nanstd(MUAv))>5;
    MUAv(out) = NaN;
    buf = reshape(MUAv,size(buf,1),size(buf,2));
    MUA(f,:) = buf;
    
    %Remove extreme trials
    meant = nanmean(MUA(f,:),2);
    out = abs((meant-nanmean(meant))./std(meant))>4;
    MUA(f(out),:) = NaN;
    
    %Kill based on pre-trial Variance
    V = var(MUA(f,bf)');
    out = abs((V-mean(V))./std(V))>4;
    MUA(f(out),:) = NaN;
    
    %get the mean background level
    back = nanmean(nanmean(MUA(f,bf)));
    %Subtract from MUA
    MUA(f,:) = MUA(f,:)-back;
    
    %Calculate the normalisation factor
    %It is quicker to first take the mean
    %of all conditions, then smooth and calculate the maximum then to
    %smooth all trials individually, this shoudl give the same result.
    %Also only take the correct trials
    f2 = find(hit(f) == 1);
    buf = nanmean(MUA(f(f2),:));
    bufs = smooth(buf,10);
    normfactor = max(bufs(pf));
    
    %Now divide by the max of the mean response
    MUA(f,:) = MUA(f,:)./normfactor;
    
    %get rid of shit days
    c = max(MUA(f,:)')-min(MUA(f,:)');
    stdcheck = nanstd(c);
    if stdcheck > 50;%2    7 6
%         keyboard;
        eraseday(FL,chan) = 1;
        MUA(f,:) = NaN(length(f),size(MUA(f,:),2));
        disp(['Warning! Erasing day ' num2str(FL) ' from channel ',num2str(chan)])
    end
    
end
CleanData = MUA';
