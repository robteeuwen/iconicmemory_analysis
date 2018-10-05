function [extractdir, sessions, rawdir, SaveDir] = infoDir_DarwinIcoMemo;

%directory of extracted data and stim log data
extractdir = 'D:\Projects\IconicMemory\datasets\Darwin\extractdata\';
% 1st relevant session : 20170327 (no 31)
sessions = {'20160606', '20161114','20170119', '20170120', '20170123', '20170124', '20170125', '20170126', '20170208', '20170209', ...
    '20170214', '20170216', '20170217', '20170220', '20170221', '20170222', '20170223', '20170224', '20170301', '20170302', '20170303',...
    '20170306', '20170308', '20170309', '20170314', '20170315', '20170316', '20170317', '20170322', '20170324',...
    '20170327', '20170328', '20170329', '20170403', '20170404', '20170410', '20170411', '20170412', '20170413', '20170418',...
    '20170419', '20170420', '20170421', '20170424', '20170425', '20170426', '20170501', '20170502', '20170503', '20170504',...
    '20170509', '20170510', '20170511', '20170512', '20170515', '20170516', '20170517', '20170519', '20170522',...
    '20170809', '20170810', '20170811', '20170816', '20170817', '20170821', '20170822', '20170823', '20170829',...
    '20170830', '20170904', '20170905', '20170906', '20170918', '20170919', '20170920', '20170921', '20171012',...
    '20171013', '20171020', '20171023', '20171024', '20171025', '20171101', '20171106', '20171108', '20171109',...
    '20171110', '20171114', '20171115', '20171116', '20171117', '20171121', '20171123', '20171124', '20171127',...
    '20180516', '20180517', '20180518', '20180524', '20180525', '20180528', '20180529', '20180530', '20180531',...
    '20180601', '20180604', '20180606', '20180607', '20180613', '20180614', '20180620', '20180621', '20180622',...
    '20180625', '20180626', '20180627', '20180709', '20180710', '20180711', '20180712', '20180717', '20180718',...
    '20180719', '20180720', '20180723', '20180724', '20180725', '20180726', '20180730', '20180731', '20180808',...
    '20180809', '20180810', '20180813', '20180814', '20180815', '20180816', '20180817', '20180820', '20180821',...
    '20180822', '20180823', '20180827', '20180828', '20180829', '20180830', '20180831', '20180903', '20180904',...
    '20180905', '20180912', '20180913', '20180914', '20180924', '20180925', '20180926', '20180927', '20181002',...
    '20181003', '20181004'};
rawdir = 'D:\Projects\IconicMemory\datasets\Darwin\logs\';
SaveDir = 'D:\Projects\IconicMemory\datasets\Darwin\Figures\';