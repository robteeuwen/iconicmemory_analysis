



rawdir = 'D:\Projects\IconicMemory\datasets\Darwin\RawData\Darwin_20180910\';

blocks = [1 3 5 7];

for block = blocks
    b = ['block-' num2str(block)];
    
    clear EVENT
    EVENT.Mytank = [rawdir];
    EVENT.Myblock = b;
    EVENT = Exinf4_all(EVENT);
    
    EVENT.Myevent = 'ENV1';
    EVENT.CHAN = 25;
    Env = Exd_all(EVENT);
    
end



