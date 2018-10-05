function [CleanData, eraseday, isgood] = NormalizeActiv(MUA, chan, days,hit ,bf, pf);



% normalize per day based on average response to connected trials
CleanData =[];
eraseday = zeros(1, numel(unique(days)));
isgood = zeros(1, max(days));
for d = unique(days)
    d_index = find(days==d);
    MUA_d = MUA(:,d_index);
    buf = MUA_d;
    MUAv = reshape(buf,size(buf,1)*size(buf,2),1);
    out = abs((MUAv-nanmean(MUAv))./nanstd(MUAv))>5;
    MUAv(out) = NaN;
    out = abs((MUAv-nanmean(MUAv))./nanstd(MUAv))>5;
    MUAv(out) = NaN;
    buf = reshape(MUAv,size(buf,1),size(buf,2));
    MUA_d = buf;
    
    % for each trial, substract baseline activity
    BL = nanmean(MUA_d(bf,:));
    MUA_d= MUA_d- ones(size(MUA_d,1),1)*BL;
    
    % remove trials if crazy variance in bl
    BLvar = nanstd(MUA_d(bf,:));
    m_var = nanstd(BLvar);
    MUA_d(:,(BLvar-nanmean(BLvar))>5*m_var) = NaN;
    
    
    
    
    
    
    hit_d = hit(d_index);
    Peak_activ = nanmean(MUA_d(pf,hit_d==1),2);
    Norm_factor = max(Peak_activ);
    MUA_d = MUA_d/Norm_factor;
    
    
    %get rid of shit days
    c = nanmax(MUA_d)-nanmin(MUA_d);
    stdcheck = nanstd(c);
    if stdcheck > 20;%2    7 6
%         keyboard;
        eraseday(unique(days)==d) = 1;
        MUA_d = NaN*MUA_d;
        disp(['Warning! Erasing day ' num2str(d) ' from channel ',num2str(chan)])
    end
    SNR = 1/nanmean(nanstd(MUA_d(bf,:),0,2));
    if SNR<0.6
        eraseday(unique(days)==d) = 1;
        MUA_d = NaN*MUA_d;
        disp(['Low SNR! Erasing day ' num2str(d) ' from channel ',num2str(chan)])
    end
    if SNR>0.6
        isgood(d) = 1;
        
    end
        
    
    if isempty(CleanData)
        CleanData = MUA_d;
    else
        CleanData = [CleanData MUA_d];
    end
end




% remove trials with crazy values
m_trial = nanmax(abs(CleanData));
CleanData(:,m_trial>50) = NaN;
m_trial = nanmax(CleanData);
v_mtrial = nanstd(m_trial);
CleanData(:,m_trial>(nanmean(m_trial)+3*v_mtrial)) = NaN;


