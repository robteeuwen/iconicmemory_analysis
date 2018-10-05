function XtractData_IcoMemo

%Xtract FigGnd data from the raw files to .../XtractedData
%Make one file per channel for fast loading
%Needs to extract around the STIMULUS bit

%This can be used with rempapped or non-remapped projects
%It adapts to teh number of channels in ENV1 and ENV2
rawdir = '\\vcnin\iconicmem\datasets\Monty\RawData\';%'\\NIN477\IconicMemory\';%
addpath(rawdir)
% rawdir = 'Z:\Daan\Bobo_Gestalt\';
xtractdir = '\\vcnin\iconicmem\datasets\Monty\extractdata\20180424\';
mkdir(xtractdir)

% info = Duvel_Logfile_Iconic;
conn2 = 0;

spikeon = 0;
lfpon = 0;
eyeon = 1;
emp = 1;
for N = 1
    
    Tankname = 'Monty20180424';
    Blockno = [1];
    
    %make a filename for the xtracted data and event file
    savename = ['Xtract_',Tankname];
    saveeyename = ['XtractEye_',Tankname];
    eventname = ['EVT_',Tankname];
    
    for B = 1:length(Blockno)
        
        blocknames = ['block-',num2str(Blockno(B))];
        
        clear EVENT
        EVENT.Mytank = [rawdir,Tankname];
        EVENT.Myblock = blocknames;
        EVENT = Exinf4_stimnew(EVENT);
        
        EVENT.Triallngth =  1.8;%1.5
        EVENT.Start =      -0.5;%-0.2
        EVENT.type = 'strms';
        
        Valid_trials = find(EVENT.Trials.target_onset>0);
%         for trial_type = 1:4
%             for corr = 1:2
%             
        %Get stimulus onset times for extraction (i.e. onset of texture)
        Trials = EVENT.Trials.stim_onset(Valid_trials);
        
        save([xtractdir,eventname,'_Block-',num2str(Blockno(B))],'EVENT')
        
        %how many channels in ENV1 and ENV2?
        for e = 1:length(EVENT.strms)
            fld1(e) = ~isempty(strmatch(EVENT.strms(e).name,'ENV1'));
%             fld2(e) = ~isempty(strmatch(EVENT.strms(e).name,'ENV2'));
        end
        e1chn = EVENT.strms(find(fld1)).channels;
%         e2chn = EVENT.strms(find(fld2)).channels;
        
        if spikeon
            clear fld1,clear fld2,clear fld3,clear fld4
            for e = 1:length(EVENT.snips)
                fld1(e) = ~isempty(strmatch(EVENT.snips(e).name,'SNP1'));
                fld2(e) = ~isempty(strmatch(EVENT.snips(e).name,'SNP2'));
                fld3(e) = ~isempty(strmatch(EVENT.snips(e).name,'SNP3'));
                %             fld4(e) = ~isempty(strmatch(EVENT.snips(e).name,'SNP4'));
            end
            s1chn = EVENT.snips(find(fld1)).channels;
            s2chn = EVENT.snips(find(fld2)).channels;
            s3chn = EVENT.snips(find(fld3)).channels;
        end
        
        
        %Save individual channels for quick loding
        for ch = 1:48%(e1chn)%+e2chn
            if ch <=e1chn
                EVENT.Myevent = 'ENV1';
                EVENT.CHAN = ch;
                Env = Exd4(EVENT, Trials);
            else
                EVENT.Myevent = 'ENV2';
                EVENT.CHAN = ch-e1chn;
                Env = Exd4(EVENT, Trials);
            end      
            if lfpon
                if ch <=e1chn
                    EVENT.Myevent = 'LFP1';
                    EVENT.CHAN = ch;
                    LFP = Exd4(EVENT, Trials);
                else
                    EVENT.Myevent = 'LFP2';
                    EVENT.CHAN = ch-e1chn;
                    LFP = Exd4(EVENT, Trials);
                end
            end
            
            if spikeon
                %SNIPS
                %Hard-coded as SNIPS can be missing if there are no snip channels
                EVENT.type = 'snips';
                if ch <=38
                    EVENT.Myevent = 'SNP1';
                    EVENT.CHAN = ch;
                    Snip = Exsnip1(EVENT, Trials);
                elseif ch <= 76
                    EVENT.Myevent = 'SNP2';
                    EVENT.CHAN = ch-38;
                    Snip = Exsnip1(EVENT, Trials);
                else
                    EVENT.Myevent = 'SNP3';
                    EVENT.CHAN = ch-76;
                    Snip = Exsnip1(EVENT, Trials);
                end
                
                %Check to see if Snip is empty
                clear j
                for i = 1:length(Snip)
                    j(i) = length(Snip{i});
                end
                emp = sum(j)==0;
            end
            
            %Save out
            if emp
                if ~lfpon
                    save([xtractdir,savename,'_Block-',num2str(Blockno(B)),'_',num2str(ch)],'Env')
                else
                    save([xtractdir,savename,'_Block-',num2str(Blockno(B)),'_',num2str(ch)],'Env','LFP')
                end
            else
                %Just save out spiketimes as otherwise files can get
                %massive!
                for g = 1:length(Snip)
                    %Extract all the snip times on this trial to a
                    %structure
                    buf = struct2cell(Snip{g}(:));
                    SnipT(g).Times =cell2mat(buf(1,:));
                end
                
                disp(['Channel ',num2str(ch),' has spikes'])
                if ~lfpon
                    save([xtractdir,savename,'_Block-',num2str(Blockno(B)),'_',num2str(ch)],'Env','SnipT')
                else
                    save([xtractdir,savename,'_Block-',num2str(Blockno(B)),'_',num2str(ch)],'Env','LFP','SnipT')
                end
            end
            disp([Tankname,'_',blocknames,'_',num2str(ch)])
        end
        %EYES?
        if eyeon
            for ch = 1:4
                EVENT.Myevent = 'Eye_';
                EVENT.CHAN = ch;
                Eye = Exd4(EVENT, Trials);
                save([xtractdir,saveeyename,'_Block-',num2str(Blockno(B)),'_',num2str(ch)],'Eye')
            end
        end
    end
end



