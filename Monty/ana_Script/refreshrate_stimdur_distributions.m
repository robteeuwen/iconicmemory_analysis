

stimdurs = unique(LOG.STIMDUR)


figure

for s = 1:length(stimdurs)
    cs = stimdurs(s);
    ct = find(LOG.STIMDUR==cs);
    
    actdur = LOG.stim_off(ct)-LOG.stim_on(ct);
    
    subplot(2,3,s)
    h{s} = histogram(actdur,30);
    title(['requested stim dur: ' num2str(cs)]);
end


maskonsets = unique(LOG.MASK_TIME);

figure
clear h
for s = 1:length(maskonsets)
    cs = maskonsets(s);
    ct = find(LOG.MASK_TIME==cs);
    
    actdur = LOG.mask_on(ct)-LOG.stim_off(ct);
    
    subplot(2,3,s)
    h{s} = histogram(actdur,30);
    title(['requested mask time: ' num2str(cs)]);
end


