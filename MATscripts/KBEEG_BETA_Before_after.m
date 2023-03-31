%fieldtrip
 clear;
addpath(genpath('/home/b1044271/Toolboxes/obob_ownft'));
obob_init_ft;

%EEGLab
addpath(genpath('/home/b1044271/Toolboxes/eeglab14_1_1b'));
eeglab; close gcf;


%% load files
cd '/home/b1044271/EEG_KB/AG_REG_Raw'
path_R='/home/b1044271/EEG_KB/AG_REG_Raw/';
dataRow_r   =struct2cell(dir(fullfile(path_R , '*.cnt'))');
files_R =dataRow_r(1,:);
% savepath= 'C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\AG_REG_Raw\preprocessed_normal_reref\';

cd '/home/b1044271/EEG_KB/AG_INV_Raw'
path_I='/home/b1044271/EEG_KB/AG_INV_Raw/';
dataRow_i   =struct2cell(dir(fullfile(path_I , '*.cnt'))');
files_I =dataRow_i(1,:);

%% Groups 
AG = [1:10, 22, 23,24,25,31,32,33];
SG = [11:21, 26,27,28,29,30];

EEG1_d_pre = nan(200,33);EEG1_d_pos = nan(200,33);
EEG2_d_pre = nan(200,33);EEG2_d_pos = nan(200,33);
EEG3_d_pre = nan(200,33);EEG3_d_pos = nan(200,33);
EEG4_d_pre = nan(200,33);EEG4_d_pos = nan(200,33);
EEG5_d_pre = nan(200,33);EEG5_d_pos = nan(200,33);
EEG6_d_pre = nan(200,33);EEG6_d_pos = nan(200,33);
EEG7_d_pre = nan(200,33);EEG7_d_pos = nan(200,33);

%% First and last
for i = 1:33

    if i <10%read data (replace with arnod delorms new function)
        EEG1=pop_loadcnt([path_R sprintf('VP0%dREG1.cnt', i)]);
        EEG2=pop_loadcnt([path_R sprintf('VP0%dREG2.cnt', i)]);
        EEG3=pop_loadcnt([path_I sprintf('VP0%dINV1.cnt', i)]);
        EEG4=pop_loadcnt([path_I sprintf('VP0%dINV2.cnt', i)]);
        EEG5=pop_loadcnt([path_I sprintf('VP0%dINV3.cnt', i)]);
        EEG6=pop_loadcnt([path_I sprintf('VP0%dINV4.cnt', i)]);
        EEG7=pop_loadcnt([path_I sprintf('VP0%dINV5.cnt', i)]);
    else
        EEG1=pop_loadcnt([path_R sprintf('VP%dREG1.cnt', i)]);
        EEG2=pop_loadcnt([path_R sprintf('VP%dREG2.cnt', i)]);
        EEG3=pop_loadcnt([path_I sprintf('VP%dINV1.cnt', i)]);
        EEG4=pop_loadcnt([path_I sprintf('VP%dINV2.cnt', i)]);
        EEG5=pop_loadcnt([path_I sprintf('VP%dINV3.cnt', i)]);
        EEG6=pop_loadcnt([path_I sprintf('VP%dINV4.cnt', i)]);
        EEG7=pop_loadcnt([path_I sprintf('VP%dINV5.cnt', i)]);
    end
    
    %% Prepro
        EEG1 = pop_select(EEG1, 'channel', 1:10); EEG2 = pop_select(EEG2, 'channel', 1:10); 
    EEG3 = pop_select(EEG3, 'channel', 1:10); EEG4 = pop_select(EEG4, 'channel', 1:10); 
    EEG5 = pop_select(EEG5, 'channel', 1:10); EEG6 = pop_select(EEG6, 'channel', 1:10); 
    EEG7 = pop_select(EEG7, 'channel', 1:10); 
      
    %add the REF channel back and perform the average reference  
     %add the REF channel back and perform the average reference   
     EEG1.nbchan                            = EEG1.nbchan+1;
    EEG1.data(end+1,:)                     = zeros(1, EEG1.pnts);
    EEG1.chanlocs(1,EEG1.nbchan).labels    = 'CZ';
    EEG1                                   = pop_chanedit(EEG1, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG1                                   = eeg_checkset( EEG1);
    EEG1                                   = pop_reref( EEG1, []);  % normal eeglab reref    
    EEG1                                   = pop_eegfiltnew(EEG1, 0.1, 60); %pop_eegfiltnew(EEG1, 13, 30);
    EEG1 = pop_cleanline(EEG1, 'bandwidth',2,'chanlist',[1:11] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    
    EEG2.nbchan                            = EEG2.nbchan+1;
    EEG2.data(end+1,:)                     = zeros(1, EEG2.pnts);
    EEG2.chanlocs(1,EEG2.nbchan).labels    = 'CZ';
    EEG2                                   = pop_chanedit(EEG1, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG2                                   = eeg_checkset( EEG2);
    EEG2                                   = pop_reref( EEG2, []);  % normal eeglab reref    
    EEG2                                   = pop_eegfiltnew(EEG2, 0.1, 60);
    EEG2 = pop_cleanline(EEG2, 'bandwidth',2,'chanlist',[1:11] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

    
%     EEG3.nbchan                          = EEG3.nbchan+1;
%     EEG3.data(end+1,:)                   = zeros(1, EEG3.pnts);
%     EEG3.chanlocs(1,EEG3.nbchan).labels  = 'CZ';
%     EEG3                                 = pop_chanedit(EEG3, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
%     EEG3                                 = eeg_checkset( EEG3);
%     EEG3                                 = pop_reref( EEG3, []);  % normal eeglab reref    
%     EEG3                                 = pop_eegfiltnew(EEG3, 0.1, 60);
%     EEG3 = pop_cleanline(EEG3, 'bandwidth',2,'chanlist',[1:11] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
% 
%     
%     EEG4.nbchan                          = EEG4.nbchan+1;
%     EEG4.data(end+1,:)                   = zeros(1, EEG4.pnts);
%     EEG4.chanlocs(1,EEG4.nbchan).labels  = 'CZ';
%     EEG4                                 = pop_chanedit(EEG4, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
%     EEG4                                 = eeg_checkset( EEG4);
%     EEG4                                 = pop_reref( EEG4, []);  % normal eeglab reref    
%     EEG4                                 = pop_eegfiltnew(EEG4, 0.1, 60);
%     EEG4 = pop_cleanline(EEG4, 'bandwidth',2,'chanlist',[1:11] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

    
    EEG5.nbchan                          = EEG5.nbchan+1;
    EEG5.data(end+1,:)                   = zeros(1, EEG5.pnts);
    EEG5.chanlocs(1,EEG5.nbchan).labels  = 'CZ';
    EEG5                                 = pop_chanedit(EEG5, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG5                                 = eeg_checkset( EEG5);
    EEG5                                 = pop_reref( EEG5, []);  % normal eeglab reref    
    EEG5                                 = pop_eegfiltnew(EEG5, 0.1, 60);
    EEG5 = pop_cleanline(EEG5, 'bandwidth',2,'chanlist',[1:11] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

    
    EEG6.nbchan                          = EEG6.nbchan+1;
    EEG6.data(end+1,:)                   = zeros(1, EEG6.pnts);
    EEG6.chanlocs(1,EEG6.nbchan).labels  = 'CZ';
    EEG6                                 = pop_chanedit(EEG6, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG6                                 = eeg_checkset( EEG6);
    EEG6                                 = pop_reref( EEG6, []);  % normal eeglab reref    
    EEG6                                 = pop_eegfiltnew(EEG6, 0.1, 60);
    EEG6 = pop_cleanline(EEG6, 'bandwidth',2,'chanlist',[1:11] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

    
    EEG7.nbchan                          = EEG7.nbchan+1;
    EEG7.data(end+1,:)                   = zeros(1, EEG7.pnts);
    EEG7.chanlocs(1,EEG7.nbchan).labels  = 'CZ';
    EEG7                                 = pop_chanedit(EEG7, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG7                                 = eeg_checkset( EEG7);
    EEG7                                 = pop_reref( EEG7, []);  % normal eeglab reref    
    EEG7                                 = pop_eegfiltnew(EEG7, 0.1, 60);
    EEG7 = pop_cleanline(EEG7, 'bandwidth',2,'chanlist',[1:11] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
   
    %% epoching - PRE (zero is the first typed letter of the word)
    
    X = find([EEG1.event(:).type] ==11);
    Xb = X(1:end-1)+2;
    preEEG1 = pop_epoch(EEG1, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG2.event(:).type] ==11);
    Xb = X(1:end-1)+2;
    preEEG2 = pop_epoch(EEG2, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG3.event(:).type] ==11);
    Xb = X(1:end-1)+2;
    preEEG3 = pop_epoch(EEG3, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
        X = find([EEG4.event(:).type] ==11);
    Xb = X(1:end-1)+2;
    preEEG4 = pop_epoch(EEG4, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG5.event(:).type] ==11);
    Xb = X(1:end-1)+2;
    preEEG5 = pop_epoch(EEG5, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG6.event(:).type] ==11);
    Xb = X(1:end-1)+2;
    preEEG6 = pop_epoch(EEG6, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG7.event(:).type] ==11);
    Xb = X(1:end-1)+2;
    preEEG7 = pop_epoch(EEG7, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    %% Epoching Post (zero is the last letter typed)
    
    X = find([EEG1.event(:).type] ==11);
    Xb = X(2:end-1)-1;
    postEEG1 = pop_epoch(EEG1, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG2.event(:).type] ==11);
    Xb = X(2:end-1)-1;
    postEEG2 = pop_epoch(EEG2, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG3.event(:).type] ==11);
    Xb = X(2:end-1)-1;
    postEEG3 = pop_epoch(EEG3, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG4.event(:).type] ==11);
    Xb = X(2:end-1)-1;
    postEEG4 = pop_epoch(EEG4, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];

    X = find([EEG5.event(:).type] ==11);
    Xb = X(2:end-1)-1;
    postEEG5 = pop_epoch(EEG5, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG6.event(:).type] ==11);
    Xb = X(2:end-1)-1;
    postEEG6 = pop_epoch(EEG6, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
    
    X = find([EEG7.event(:).type] ==11);
    Xb = X(2:end-1)-1;
    postEEG7 = pop_epoch(EEG7, {}, [-5 5], 'eventindices',Xb);
    X=[]; Xb=[];
   
    %% Find the duration of 1)fixation cross,  
    X  = find([EEG5.event(:).type] ==11); X=X(1:end-1);
    X1 = [EEG5.event(X).latency]; X2 = [EEG5.event(X+1).latency];
    X3 = (X2-X1)./500; % duration of fixation corss
    X=[];X1=[];X2=[];X3=[];
    
    %% Duration of 2) ISI between word-presentation and typing first letter
        X  = find([EEG1.event(:).type] ==11); X=X(1:end-1);
        X1 = [EEG1.event(X+1).latency]; X2 = [EEG1.event(X+2).latency];
        EEG1_d_pre(1:length(X1),i) = (X2-X1)./500; X=[];X1=[];X2=[];
        
        X  = find([EEG2.event(:).type] ==11); X=X(1:end-1);
        X1 = [EEG2.event(X+1).latency]; X2 = [EEG2.event(X+2).latency];
        EEG2_d_pre(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];
        
        X  = find([EEG3.event(:).type] ==11); X=X(1:end-1);
        X1 = [EEG3.event(X+1).latency]; X2 = [EEG3.event(X+2).latency];
        EEG3_d_pre(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];
        
        X  = find([EEG4.event(:).type] ==11); X=X(1:end-1);
        X1 = [EEG4.event(X+1).latency]; X2 = [EEG4.event(X+2).latency];
        EEG4_d_pre(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];

        X  = find([EEG5.event(:).type] ==11); X=X(1:end-1);
        X1 = [EEG5.event(X+1).latency]; X2 = [EEG5.event(X+2).latency];
        EEG5_d_pre(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];
                        
        X  = find([EEG6.event(:).type] ==11); X=X(1:end-1);
        X1 = [EEG6.event(X+1).latency]; X2 = [EEG6.event(X+2).latency];
        EEG6_d_pre(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];               
        
        X  = find([EEG7.event(:).type] ==11); X=X(1:end-1);
        X1 = [EEG7.event(X+1).latency]; X2 = [EEG7.event(X+2).latency];
        EEG7_d_pre(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];
    
    %% 3) Duration between last letter and 1st letter of the following word (helps for beta rebound calc)
        X  = find([EEG1.event(:).type] ==11); X=X(2:end-1);
        X1 = [EEG1.event(X-1).latency]; X2 = [EEG1.event(X+2).latency];
        EEG1_d_pos(1:length(X1),i) = (X2-X1)./500; X=[];X1=[];X2=[];
        
        X  = find([EEG2.event(:).type] ==11); X=X(2:end-1);
        X1 = [EEG2.event(X-1).latency]; X2 = [EEG2.event(X+2).latency];
        EEG2_d_pos(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];
        
        X  = find([EEG3.event(:).type] ==11); X=X(2:end-1);
        X1 = [EEG3.event(X-1).latency]; X2 = [EEG3.event(X+2).latency];
        EEG3_d_pos(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];
        
        X  = find([EEG4.event(:).type] ==11); X=X(2:end-1);
        X1 = [EEG4.event(X-1).latency]; X2 = [EEG4.event(X+2).latency];
        EEG4_d_pos(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];

        X  = find([EEG5.event(:).type] ==11); X=X(2:end-1);
        X1 = [EEG5.event(X-1).latency]; X2 = [EEG5.event(X+2).latency];
        EEG5_d_pos(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];
                        
        X  = find([EEG6.event(:).type] ==11); X=X(2:end-1);
        X1 = [EEG6.event(X-1).latency]; X2 = [EEG6.event(X+2).latency];
        EEG6_d_pos(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];               
        
        X  = find([EEG7.event(:).type] ==11); X=X(2:end-1);
        X1 = [EEG7.event(X-1).latency]; X2 = [EEG7.event(X+2).latency];
        EEG7_d_pos(1:length(X1),i) = (X2-X1)./500;X=[];X1=[];X2=[];
    
    %% run bad epoch detection and reject bad epoches
    [preEEG1, rmepochs] = pop_autorej( preEEG1,'nogui','on');epoch_props = epoch_properties(preEEG1,1:preEEG1.nbchan);
    BadEpochs   = min_z(epoch_props);preEEG1 = pop_rejepoch(preEEG1, BadEpochs,0);
    
    [preEEG2, rmepochs] = pop_autorej(preEEG2,'nogui','on');epoch_props = epoch_properties(preEEG2,1:preEEG2.nbchan);
    BadEpochs   = min_z(epoch_props);preEEG2 = pop_rejepoch(preEEG2, BadEpochs,0);

    [preEEG3, rmepochs] = pop_autorej( preEEG3,'nogui','on');epoch_props = epoch_properties(preEEG3,1:preEEG3.nbchan);
    BadEpochs   = min_z(epoch_props);preEEG3 = pop_rejepoch(preEEG3, BadEpochs,0);
    
    [preEEG4, rmepochs] = pop_autorej( preEEG4,'nogui','on');epoch_props = epoch_properties(preEEG4,1:preEEG4.nbchan);
    BadEpochs   = min_z(epoch_props);preEEG4 = pop_rejepoch(preEEG4, BadEpochs,0);
    
    [preEEG5, rmepochs] = pop_autorej( preEEG5,'nogui','on');epoch_props = epoch_properties(preEEG5,1:preEEG5.nbchan);
    BadEpochs   = min_z(epoch_props);preEEG5 = pop_rejepoch(preEEG5, BadEpochs,0);
    
    [preEEG6, rmepochs] = pop_autorej( preEEG6,'nogui','on');epoch_props = epoch_properties(preEEG6,1:preEEG6.nbchan);
    BadEpochs   = min_z(epoch_props);preEEG6 = pop_rejepoch(preEEG6, BadEpochs,0);
    
    [preEEG7, rmepochs] = pop_autorej(preEEG7,'nogui','on');epoch_props = epoch_properties(preEEG7,1:preEEG7.nbchan);
    BadEpochs   = min_z(epoch_props);preEEG7 = pop_rejepoch(preEEG7, BadEpochs,0); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% POST LAST LETTER WORD TYPED
        [postEEG1, rmepochs] = pop_autorej( postEEG1,'nogui','on');epoch_props = epoch_properties(postEEG1,1:postEEG1.nbchan);
    BadEpochs   = min_z(epoch_props);EEG1 = pop_rejepoch(postEEG1, BadEpochs,0);
    
    [postEEG2, rmepochs] = pop_autorej(postEEG2,'nogui','on');epoch_props = epoch_properties(postEEG2,1:postEEG2.nbchan);
    BadEpochs   = min_z(epoch_props);postEEG2 = pop_rejepoch(postEEG2, BadEpochs,0);

    [postEEG3, rmepochs] = pop_autorej( postEEG3,'nogui','on');epoch_props = epoch_properties(postEEG3,1:postEEG3.nbchan);
    BadEpochs   = min_z(epoch_props);postEEG3 = pop_rejepoch(postEEG3, BadEpochs,0);
    
    [postEEG4, rmepochs] = pop_autorej(postEEG4,'nogui','on');epoch_props = epoch_properties(postEEG4,1:postEEG4.nbchan);
    BadEpochs   = min_z(epoch_props);postEEG4 = pop_rejepoch(postEEG4, BadEpochs,0);
    
    [postEEG5, rmepochs] = pop_autorej(postEEG5,'nogui','on');epoch_props = epoch_properties(postEEG5,1:postEEG5.nbchan);
    BadEpochs   = min_z(epoch_props);postEEG5 = pop_rejepoch(postEEG5, BadEpochs,0);
    
    [postEEG6, rmepochs] = pop_autorej(postEEG6,'nogui','on');epoch_props = epoch_properties(postEEG6,1:postEEG6.nbchan);
    BadEpochs   = min_z(epoch_props);postEEG6 = pop_rejepoch(postEEG6, BadEpochs,0);
    
    [postEEG7, rmepochs] = pop_autorej(postEEG7,'nogui','on');epoch_props = epoch_properties(postEEG7,1:postEEG7.nbchan);
    BadEpochs   = min_z(epoch_props);postEEG7 = pop_rejepoch(postEEG7, BadEpochs,0);
    
    %% Equal number of epochs
    XR = min([length(preEEG1.epoch) length(preEEG2.epoch)]);
    XI = min([length(preEEG3.epoch) length(preEEG4.epoch) length(preEEG5.epoch) length(preEEG6.epoch)]);
    
    XR2 = min([length(postEEG1.epoch) length(postEEG2.epoch)]);
    XI2 = min([length(postEEG3.epoch) length(postEEG4.epoch) length(postEEG5.epoch) length(postEEG6.epoch)]);
    
    preEEG1a  = pop_select(preEEG1,'trial',datasample(1:length(preEEG1.epoch),XR,'Replace',false));
    preEEG2a  = pop_select(preEEG2,'trial',datasample(1:length(preEEG2.epoch),XR,'Replace',false));
    preEEG3a  = pop_select(preEEG3,'trial',datasample(1:length(preEEG3.epoch),XI,'Replace',false));
    preEEG4a  = pop_select(preEEG4,'trial',datasample(1:length(preEEG4.epoch),XI,'Replace',false));
    preEEG5a  = pop_select(preEEG5,'trial',datasample(1:length(preEEG5.epoch),XI,'Replace',false));
    preEEG6a  = pop_select(preEEG6,'trial',datasample(1:length(preEEG6.epoch),XI,'Replace',false));
    
    postEEG1a  = pop_select(postEEG1,'trial',datasample(1:length(postEEG1.epoch),XR2,'Replace',false));
    postEEG2a  = pop_select(postEEG2,'trial',datasample(1:length(postEEG2.epoch),XR2,'Replace',false));
    postEEG3a  = pop_select(postEEG3,'trial',datasample(1:length(postEEG3.epoch),XI2,'Replace',false));
    postEEG4a  = pop_select(postEEG4,'trial',datasample(1:length(postEEG4.epoch),XI2,'Replace',false));
    postEEG5a  = pop_select(postEEG5,'trial',datasample(1:length(postEEG5.epoch),XI2,'Replace',false));
    postEEG6a  = pop_select(postEEG6,'trial',datasample(1:length(postEEG6.epoch),XI2,'Replace',false));

   preEEG1b = eeglab2fieldtrip(preEEG1a,'preprocessing','none');
   preEEG2b = eeglab2fieldtrip(preEEG2a,'preprocessing','none');
   preEEG3b = eeglab2fieldtrip(preEEG3a,'preprocessing','none');
   preEEG4b = eeglab2fieldtrip(preEEG4a,'preprocessing','none');
   preEEG5b = eeglab2fieldtrip(preEEG5a,'preprocessing','none');
   preEEG6b = eeglab2fieldtrip(preEEG6a,'preprocessing','none');
   preEEG7b = eeglab2fieldtrip(preEEG7,'preprocessing','none'); 
      
   postEEG1b = eeglab2fieldtrip(postEEG1a,'preprocessing','none');
   postEEG2b = eeglab2fieldtrip(postEEG2a,'preprocessing','none');
   postEEG3b = eeglab2fieldtrip(postEEG3a,'preprocessing','none');
   postEEG4b = eeglab2fieldtrip(postEEG4a,'preprocessing','none');
   postEEG5b = eeglab2fieldtrip(postEEG5a,'preprocessing','none');
   postEEG6b = eeglab2fieldtrip(postEEG6a,'preprocessing','none');
   postEEG7b = eeglab2fieldtrip(postEEG7,'preprocessing','none');
   
    %% TF transformation
    
  cfg=[];
  cfg.method       = 'mtmconvol'; 
  cfg.taper           = 'dpss';
  cfg.output         =  'pow';
  cfg.foi                =  0.5:0.5:45;
  cfg.t_ftimwin     = 5./cfg.foi;  % 5 cycles per time window
  cfg.toi                = -5:0.01:5;
  cfg.tapsmofrq   = 0.4*cfg.foi; % MULTI-TAPER
  cfg.keeptrials     = 'yes';
  
  PREEG1{i} = ft_freqanalysis(cfg, preEEG1b);
  PREEG2{i} = ft_freqanalysis(cfg, preEEG2b);
  PREEG3{i} = ft_freqanalysis(cfg, preEEG3b);
  PREEG4{i} = ft_freqanalysis(cfg, preEEG4b);
  PREEG5{i} = ft_freqanalysis(cfg, preEEG5b);
  PREEG6{i} = ft_freqanalysis(cfg, preEEG6b);
  PREEG7{i} = ft_freqanalysis(cfg, preEEG7b);

  POEEG1{i} = ft_freqanalysis(cfg, postEEG1b);
  POEEG2{i}= ft_freqanalysis(cfg, postEEG2b);
  POEEG3{i} = ft_freqanalysis(cfg, postEEG3b);
  POEEG4{i} = ft_freqanalysis(cfg, postEEG4b);
  POEEG5{i} = ft_freqanalysis(cfg, postEEG5b);
  POEEG6{i} = ft_freqanalysis(cfg, postEEG6b);
  POEEG7{i} = ft_freqanalysis(cfg, postEEG7b);
  
    
end

cfg_g =[];
for i = 1:33

PRE1{i}=ft_freqdescriptives(cfg_g, PREEG1{i});
PRE2{i}=ft_freqdescriptives(cfg_g, PREEG2{i});
PRE3{i}=ft_freqdescriptives(cfg_g, PREEG3{i});
PRE4{i}=ft_freqdescriptives(cfg_g, PREEG4{i});
PRE5{i}=ft_freqdescriptives(cfg_g, PREEG5{i});
PRE6{i}=ft_freqdescriptives(cfg_g, PREEG6{i});
PRE7{i}=ft_freqdescriptives(cfg_g, PREEG7{i});

POS1{i}=ft_freqdescriptives(cfg_g, POEEG1{i});
POS2{i}=ft_freqdescriptives(cfg_g, POEEG2{i});
POS3{i}=ft_freqdescriptives(cfg_g, POEEG3{i});
POS4{i}=ft_freqdescriptives(cfg_g, POEEG4{i});
POS5{i}=ft_freqdescriptives(cfg_g, POEEG5{i});
POS6{i}=ft_freqdescriptives(cfg_g, POEEG6{i});
POS7{i}=ft_freqdescriptives(cfg_g, POEEG7{i});

end

%% PRE
cfg_g=[];
cfg_g.keepindividual = 'yes';

PRE1_S=ft_freqgrandaverage(cfg_g, PRE1{SG});
PRE2_S=ft_freqgrandaverage(cfg_g, PRE2{SG});
PRE3_S=ft_freqgrandaverage(cfg_g, PRE3{SG});
PRE4_S=ft_freqgrandaverage(cfg_g, PRE4{SG});
PRE5_S=ft_freqgrandaverage(cfg_g, PRE5{SG});
PRE6_S=ft_freqgrandaverage(cfg_g, PRE6{SG});
PRE7_S=ft_freqgrandaverage(cfg_g, PRE7{SG});

PRE1_W=ft_freqgrandaverage(cfg_g, PRE1{AG});
PRE2_W=ft_freqgrandaverage(cfg_g, PRE2{AG});
PRE3_W=ft_freqgrandaverage(cfg_g, PRE3{AG});
PRE4_W=ft_freqgrandaverage(cfg_g, PRE4{AG});
PRE5_W=ft_freqgrandaverage(cfg_g, PRE5{AG});
PRE6_W=ft_freqgrandaverage(cfg_g, PRE6{AG});
PRE7_W=ft_freqgrandaverage(cfg_g, PRE7{AG});




%% POST
PO1_S=ft_freqgrandaverage(cfg_g, POS1{SG});
PO2_S=ft_freqgrandaverage(cfg_g, POS2{SG});
PO3_S=ft_freqgrandaverage(cfg_g, POS3{SG});
PO4_S=ft_freqgrandaverage(cfg_g, POS4{SG});
PO5_S=ft_freqgrandaverage(cfg_g, POS5{SG});
PO6_S=ft_freqgrandaverage(cfg_g, POS6{SG});
PO7_S=ft_freqgrandaverage(cfg_g, POS7{SG});

PO1_W=ft_freqgrandaverage(cfg_g, POS1{AG});
PO2_W=ft_freqgrandaverage(cfg_g, POS2{AG});
PO3_W=ft_freqgrandaverage(cfg_g, POS3{AG});
PO4_W=ft_freqgrandaverage(cfg_g, POS4{AG});
PO5_W=ft_freqgrandaverage(cfg_g, POS5{AG});
PO6_W=ft_freqgrandaverage(cfg_g, POS6{AG});
PO7_W=ft_freqgrandaverage(cfg_g, POS7{AG});


%% STATS

cfg=[];
cfg.statistic        ='ft_statfun_depsamplesT'; %between groups independant sampels statistics
cfg.method           ='montecarlo'; % Monte Carlo method for calculating the signif prob an estimate of the p-value under the permutation distribution.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.avgoverchan      = 'yes';
% cfg.frequency        = [13 30];
% cfg.avgoverfreq      = 'yes';
% cfg.channel      = {'C3','C4'};
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
% cfg.minnbchan        = 2;
cfg.uvar             = 1;
cfg.ivar             = 2;
% Design Matrix for T-Test (2 Conditions)
subj = 16;
design = zeros(2,2*subj);
for m = 1:subj
        design(1,m) = m;
end
for m = 1:subj
        design(1,subj+m) = m;
        
end
design(2,1:subj) = 1;
design(2,subj+1:2*subj) = 2;
cfg.design= design; 

%  load('/home/b1044271/Toolboxes/LAYOUT/ChanInfoFT.mat')
% cfg.neighbours    = neighb;



cfg.latency          = [-1 3];
[st_f]  =   ft_freqstatistics(cfg, PO6_S, PO5_S);
[st_fa]  =   ft_freqstatistics(cfg, PO5_S, PO1_S);
[st_fb]  =   ft_freqstatistics(cfg, PO6_S, PO2_S);
[st_fr1]  =   ft_freqstatistics(cfg, PO2_S, PO1_S);


cfg.latency          = [-3 1];
[st_f2]  =   ft_freqstatistics(cfg, PRE6_S, PRE5_S);
[st_f2a]  =   ft_freqstatistics(cfg, PRE1_S, PRE5_S);
[st_f2b]  =   ft_freqstatistics(cfg, PRE2_S, PRE6_S);
[st_fr2]  =   ft_freqstatistics(cfg, PRE2_S, PRE1_S);


subj = 17;
design = zeros(2,2*subj);
for m = 1:subj
        design(1,m) = m;
end

for m = 1:subj
        design(1,subj+m) = m;
        
end
design(2,1:subj) = 1;
design(2,subj+1:2*subj) = 2;
cfg.design= design;       

cfg.latency          = [-1 3];
[st_fw]    =   ft_freqstatistics(cfg, PO6_W, PO5_W);
[st_fwa]  =   ft_freqstatistics(cfg, PO5_W, PO1_W);
[st_fwb]  =   ft_freqstatistics(cfg, PO6_W, PO2_W);

cfg.latency          = [-3 1];
[st_fw2]    =   ft_freqstatistics(cfg, PRE6_W, PRE5_W);
[st_fw2a]  =   ft_freqstatistics(cfg, PRE1_W, PRE5_W);
[st_fw2b]  =   ft_freqstatistics(cfg, PRE2_W, PRE6_W);

%% TOPOs

cfg = [];
cfg.zlim = [-40 40];
cfg.highlight = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = X;
cfg.colorbar= 'yes';
cfg.colormap=colormap(brewermap(256, '*RdYlBu'));   
cfg.style = 'straight_imsat';
cfg.gridscale =800;
cfg.markersize =0.5;
cfg.highlightsize=18;
cfg.highlightcolor=[0 0 0];
cfg.highlightsymbol='.';
cfg.highlightchannel = find(st_f.mask);

ft_topoplotER(cfg, st_f); colorbar; title('baseline');
[st_fb]  =   ft_freqstatistics(cfg, PO2_S, PO6_S);


%% PLOT POST-MOVEMENT Power but first average over BBA
cfg             = [];
cfg.frequency   = [13 30];
cfg.avgoverfreq = 'yes';
%  cfg.channel     = {'C3','C4'};
% cfg.avgoverchan = 'yes';
cfg.latency = [0 2];
PO1_SB =ft_selectdata(cfg, PO1_S);
PO2_SB =ft_selectdata(cfg, PO2_S);
PO3_SB =ft_selectdata(cfg, PO3_S);
PO4_SB =ft_selectdata(cfg, PO4_S);
PO5_SB =ft_selectdata(cfg, PO5_S);
PO6_SB =ft_selectdata(cfg, PO6_S);
PO7_SB =ft_selectdata(cfg, PO7_S);
%%%
PO1_WB =ft_selectdata(cfg, PO1_W);
PO2_WB =ft_selectdata(cfg, PO2_W);
PO3_WB =ft_selectdata(cfg, PO3_W);
PO4_WB =ft_selectdata(cfg, PO4_W);
PO5_WB =ft_selectdata(cfg, PO5_W);
PO6_WB =ft_selectdata(cfg, PO6_W);
PO7_WB =ft_selectdata(cfg, PO7_W);

po1s_a = squeeze(mean(PO1_SB.powspctrm));
po2s_a = squeeze(mean(PO2_SB.powspctrm));
po5s_a = squeeze(mean(PO5_SB.powspctrm));
po6s_a = squeeze(mean(PO6_SB.powspctrm));

po1s_s = squeeze(std(PO1_SB.powspctrm))./sqrt(16);
po2s_s = squeeze(std(PO2_SB.powspctrm))./sqrt(16);
po5s_s = squeeze(std(PO5_SB.powspctrm))./sqrt(16);
po6s_s = squeeze(std(PO6_SB.powspctrm))./sqrt(16);

po1w_a = squeeze(mean(PO1_WB.powspctrm));
po2w_a = squeeze(mean(PO2_WB.powspctrm));
po5w_a = squeeze(mean(PO5_WB.powspctrm));
po6w_a = squeeze(mean(PO6_WB.powspctrm));

po1w_s = squeeze(std(PO1_WB.powspctrm))./sqrt(17);
po2w_s = squeeze(std(PO2_WB.powspctrm))./sqrt(17);
po5w_s = squeeze(std(PO5_WB.powspctrm))./sqrt(17);
po6w_s = squeeze(std(PO6_WB.powspctrm))./sqrt(17);

f=figure;
colors1 = [0 0.5 0; 0.25 0.25 0.25];
subplot(1,2,1); bl(:,1)=boundedline(PO5_SB.time, po1s_a,po1s_s, PO5_SB.time, po5s_a,po5s_s, 'cmap', colors1, 'alpha','transparency', 0.2); 
set(gca,'FontName','sans serif','FontSize',20);
line([0 0], [-0.5 3.], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');

subplot(1,2,2); bl(:,2)=boundedline(PO5_WB.time, po1w_a,po1w_s, PO5_SB.time, po5w_a,po5w_s, 'cmap', colors1, 'alpha','transparency', 0.2); 
set(gca,'FontName','sans serif','FontSize',14);
line([0 0], [-0.5 3.], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');
set(bl(:),'LineStyle','-','LineWidth',4) %set line width

% ylim([0 1.5])
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
lh = legend(bl(:));
legnames = {'I3', 'I4'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})]; 
end
lh.String = str; lh.FontSize=24; lh.Box = 'off'; lh.Location = 'northwest';
ax = findobj(f,'Type','Axes');
namess = {'REG','INV'};

for i=1:length(ax)
    ylim(ax(i),[0.5  3.])
    xlim(ax(i),[-1 3])
     title(ax(i),namess{3-i})
    set(ax(i),'FontName','sans serif','FontSize',14);
%   set([ hXLabel(i) hYLabel(i) ], 'FontName','sans serif','FontSize',24);
 set(ax(i),'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
end

% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes

han=axes(f,'visible','off');
han.XLabel.Visible='on';
xlabel(han,sprintf('\n time (s)'), 'FontName','sans serif','FontSize',24);
han.YLabel.Visible='on';
ylabel(han,{'Beta Power (uv2)', ' '}, 'FontName','sans serif','FontSize',24);


%% Permutations

cfg=[];
cfg.statistic        ='ft_statfun_depsamplesT'; %between groups independant sampels statistics
cfg.method           ='montecarlo'; % Monte Carlo method for calculating the signif prob an estimate of the p-value under the permutation distribution.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.avgoverchan      = 'yes';

cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
cfg.uvar             = 1;
cfg.ivar             = 2;

% Design Matrix for T-Test (2 Conditions)
subj = 16;
design = zeros(2,2*subj);
for m = 1:subj
        design(1,m) = m;
end
for m = 1:subj
        design(1,subj+m) = m;
        
end
design(2,1:subj) = 1;
design(2,subj+1:2*subj) = 2;
cfg.design= design; 

[st_f]  =   ft_freqstatistics(cfg, PO6_SB, PO5_SB);

subj = 17;
design = zeros(2,2*subj);
for m = 1:subj
        design(1,m) = m;
end
for m = 1:subj
        design(1,subj+m) = m;
        
end
design(2,1:subj) = 1;
design(2,subj+1:2*subj) = 2;
cfg.design= design; 

[st_fw]  =   ft_freqstatistics(cfg, PO6_WB, PO5_WB);
%% TIME SERIES - PRE
cfg                        = [];
cfg.frequency       = [13 30];
cfg.avgoverfreq   = 'yes';
% cfg.channel       = {'C3','C4'};
% cfg.avgoverchan  = 'yes';
cfg.avgovertime =  'yes';
cfg.latency = [-1 0];

PRE1_SB =ft_selectdata(cfg, PRE1_S);
PRE2_SB =ft_selectdata(cfg, PRE2_S);
PRE3_SB =ft_selectdata(cfg, PRE3_S);
PRE4_SB =ft_selectdata(cfg, PRE4_S);
PRE5_SB =ft_selectdata(cfg, PRE5_S);
PRE6_SB =ft_selectdata(cfg, PRE6_S);
PRE7_SB =ft_selectdata(cfg, PRE7_S);

%pre-movement BETA per subj
writetable(table(PRE5_SB.powspctrm),'Pre-movement_Beta_SG_pre.txt', 'WriteVariableNames', 0)
writetable(table(PRE6_SB.powspctrm),'Pre-movement_Beta_SG_post.txt', 'WriteVariableNames', 0)

PRE1_WB =ft_selectdata(cfg, PRE1_W);
PRE2_WB =ft_selectdata(cfg, PRE2_W);
PRE3_WB =ft_selectdata(cfg, PRE3_W);
PRE4_WB =ft_selectdata(cfg, PRE4_W);
PRE5_WB =ft_selectdata(cfg, PRE5_W);
PRE6_WB =ft_selectdata(cfg, PRE6_W);
PRE7_WB =ft_selectdata(cfg, PRE7_W);


pe1s_a = squeeze(mean(PRE1_SB.powspctrm));
pe2s_a = squeeze(mean(PRE2_SB.powspctrm));
pe5s_a = squeeze(mean(PRE5_SB.powspctrm));
pe6s_a = squeeze(mean(PRE6_SB.powspctrm));

pe1s_s = squeeze(std(PRE1_SB.powspctrm))./sqrt(16);
pe2s_s = squeeze(std(PRE2_SB.powspctrm))./sqrt(16);
pe5s_s = squeeze(std(PRE5_SB.powspctrm))./sqrt(16);
pe6s_s = squeeze(std(PRE6_SB.powspctrm))./sqrt(16);

%%%%%%%%%%%%%%%%%%%%%

pe1w_a = squeeze(mean(PRE1_WB.powspctrm));
pe2w_a = squeeze(mean(PRE2_WB.powspctrm));
pe5w_a = squeeze(mean(PRE5_WB.powspctrm));
pe6w_a = squeeze(mean(PRE6_WB.powspctrm));

pe1w_s = squeeze(std(PRE1_WB.powspctrm))./sqrt(17);
pe2w_s = squeeze(std(PRE2_WB.powspctrm))./sqrt(17);
pe5w_s = squeeze(std(PRE5_WB.powspctrm))./sqrt(17);
pe6w_s = squeeze(std(PRE6_WB.powspctrm))./sqrt(17);

f=figure;
colors1 = [0 0.5 0; 0.25 0.25 0.25];
subplot(1,2,1); bl(:,1)=boundedline(PRE5_SB.time, pe5s_a,pe5s_s, PRE5_SB.time, pe6s_a,pe6s_s, 'cmap', colors1, 'alpha','transparency', 0.2); 
set(gca,'FontName','sans serif','FontSize',20);
line([0 0], [-0.5 3.], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');

subplot(1,2,2); bl(:,2)=boundedline(PRE5_WB.time, pe5w_a,pe5w_s, PRE5_SB.time, pe6w_a,pe6w_s, 'cmap', colors1, 'alpha','transparency', 0.2); 
set(gca,'FontName','sans serif','FontSize',14);
line([0 0], [-0.5 3.], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');


set(bl(:),'LineStyle','-','LineWidth',4) %set line width

% ylim([0 1.5])
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
lh = legend(bl(:));
legnames = {'I3', 'I4'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})]; 
end
lh.String = str; lh.FontSize=24; lh.Box = 'off'; lh.Location = 'northwest';
ax = findobj(f,'Type','Axes');
namess = {'Sleep','Wake'};

for i=1:length(ax)
    ylim(ax(i),[0.5  3])
%     xlim(ax(i),[-3 1])
     title(ax(i),namess{3-i})
    set(ax(i),'FontName','sans serif','FontSize',14);
%   set([ hXLabel(i) hYLabel(i) ], 'FontName','sans serif','FontSize',24);
 set(ax(i),'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
end

% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes

han=axes(f,'visible','off');
han.XLabel.Visible='on';
xlabel(han,sprintf('\n time (s)'), 'FontName','sans serif','FontSize',24);
han.YLabel.Visible='on';
ylabel(han,{'Beta Power (uv2)', ' '}, 'FontName','sans serif','FontSize',24);

%% POST

cfg                        = [];
cfg.frequency       = [13 30];
cfg.avgoverfreq   = 'yes';
% cfg.channel       = {'C3','C4'};
% cfg.avgoverchan  = 'yes';
cfg.latency = [0 2];
cfg.avgovertime='yes';

POS1_SB =ft_selectdata(cfg, PO1_S);
POS2_SB =ft_selectdata(cfg, PO2_S);
POS3_SB =ft_selectdata(cfg, PO3_S);
POS4_SB =ft_selectdata(cfg, PO4_S);
POS5_SB =ft_selectdata(cfg, PO5_S);
POS6_SB =ft_selectdata(cfg, PO6_S);
POS7_SB =ft_selectdata(cfg, PO7_S);

%pre-movement BETA per subj
writetable(table(POS5_SB.powspctrm),'Post-movement_Beta_SG_pre_ALL.txt', 'WriteVariableNames', 0)
writetable(table(POS6_SB.powspctrm),'Post-movement_Beta_SG_post_ALL.txt', 'WriteVariableNames', 0)

POS1_WB =ft_selectdata(cfg, PO1_W);
POS2_WB =ft_selectdata(cfg, PO2_W);
POS3_WB =ft_selectdata(cfg, PO3_W);
POS4_WB =ft_selectdata(cfg, PO4_W);
POS5_WB =ft_selectdata(cfg, PO5_W);
POS6_WB =ft_selectdata(cfg, PO6_W);
POS7_WB =ft_selectdata(cfg, PO7_W);

figure();
for i = 1:16
   subplot(5,4,i); plot(squeeze(POS1_WB.powspctrm(i,:,:,:)))
end


pe1s_a = squeeze(mean(POS1_SB.powspctrm));
pe2s_a = squeeze(mean(POS2_SB.powspctrm));
pe5s_a = squeeze(mean(POS5_SB.powspctrm));
pe6s_a = squeeze(mean(POS6_SB.powspctrm));

pe1s_s = squeeze(std(POS1_SB.powspctrm))./sqrt(16);
pe2s_s = squeeze(std(POS2_SB.powspctrm))./sqrt(16);
pe5s_s = squeeze(std(POS5_SB.powspctrm))./sqrt(16);
pe6s_s = squeeze(std(POS6_SB.powspctrm))./sqrt(16);

%%%%%%%%%%%%%%%%%%%%%

pe1w_a = squeeze(mean(POS1_WB.powspctrm));
pe2w_a = squeeze(mean(POS2_WB.powspctrm));
pe5w_a = squeeze(mean(POS5_WB.powspctrm));
pe6w_a = squeeze(mean(POS6_WB.powspctrm));

pe1w_s = squeeze(std(POS1_WB.powspctrm))./sqrt(16);
pe2w_s = squeeze(std(POS2_WB.powspctrm))./sqrt(16);
pe5w_s = squeeze(std(POS5_WB.powspctrm))./sqrt(16);
pe6w_s = squeeze(std(POS6_WB.powspctrm))./sqrt(16);

f=figure;
colors1 = [0 0.5 0; 0.25 0.25 0.25];
subplot(1,2,1); bl(:,1)=boundedline(POS5_SB.time, pe5s_a,pe5s_s, POS5_SB.time, pe1s_a,pe1s_s, 'cmap', colors1, 'alpha','transparency', 0.2); 
set(gca,'FontName','sans serif','FontSize',20);
line([0 0], [-0.5 3.], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');

subplot(1,2,2); bl(:,2)=boundedline(POS5_WB.time, pe5w_a,pe5w_s, POS5_SB.time, pe1w_a,pe1w_s, 'cmap', colors1, 'alpha','transparency', 0.2); 
set(gca,'FontName','sans serif','FontSize',14);
line([0 0], [-0.5 3.], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');


set(bl(:),'LineStyle','-','LineWidth',4) %set line width

% ylim([0 1.5])
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
lh = legend(bl(:));
legnames = {'INV3', 'REG1'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})]; 
end
lh.String = str; lh.FontSize=24; lh.Box = 'off'; lh.Location = 'northwest';
ax = findobj(f,'Type','Axes');
namess = {'Sleep','Wake'};

for i=1:length(ax)
%     ylim(ax(i),[0.5  3])
     xlim(ax(i),[-1 3])
     title(ax(i),namess{3-i})
    set(ax(i),'FontName','sans serif','FontSize',14);
%   set([ hXLabel(i) hYLabel(i) ], 'FontName','sans serif','FontSize',24);
 set(ax(i),'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
end

% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes

han=axes(f,'visible','off');
han.XLabel.Visible='on';
xlabel(han,sprintf('\n time (s)'), 'FontName','sans serif','FontSize',24);
han.YLabel.Visible='on';
ylabel(han,{'Beta Power (uv2)', ' '}, 'FontName','sans serif','FontSize',24);


%% Single trials
cfg                        = [];
cfg.frequency       = [13 30];
cfg.avgoverfreq   = 'yes';
%  cfg.channel       = {'C3', 'C4'};
% cfg.avgoverchan  = 'yes';
cfg.latency = [0 1.5];

R1m=nan(1000,33);R1a=nan(1000,33);
I3m=nan(1000,33);I3a=nan(1000,33);
R2m=nan(1000,33);R2a=nan(1000,33);
I4m=nan(1000,33);I4a=nan(1000,33);


for ii = 1:33
POX1 = ft_selectdata(cfg,PREEG1{ii});
POX2 = ft_selectdata(cfg,PREEG2{ii});
POX5 = ft_selectdata(cfg,PREEG5{ii});
POX6 = ft_selectdata(cfg,PREEG6{ii});

POX1p = squeeze(POX1.powspctrm);
for iv = 1:size(POX1p ,1)
    X=[];Xm=[];Y=[];
    X        = POX1p(iv,:);
    Xm     = mean(X)+(1.75*(std(X)));
    Y        = find(X > Xm);
   R1m(iv,ii)   = length(Y);
   R1a(iv,ii) = mean(X(Y));
end


POX2p = squeeze(POX2.powspctrm);
for iv = 1:size(POX2p ,1)
    X=[];Xm=[];Y=[];
    X        = POX2p(iv,:);
    Xm     = mean(X)+(1.75*(std(X)));
    Y        = find(X > Xm);
   R2m(iv,ii)   = length(Y);
   R2a(iv,ii) = mean(X(Y));
end

POX5p = squeeze(POX5.powspctrm);
for iv = 1:size(POX5p ,1)
    X=[];Xm=[];Y=[];
    X        = POX5p(iv,:);
    Xm     = mean(X)+(1.75*(std(X)));
    Y        = find(X > Xm);
   I3m(iv,ii)   = length(Y);
   I3a(iv,ii) = mean(X(Y));
end
I3t(ii)=length(POX5p);
POX6p = squeeze(POX6.powspctrm);
for iv = 1:size(POX6p ,1)
    X=[];Xm=[];Y=[];
    X        = POX6p(iv,:);
    Xm     = mean(X)+(1.75*(std(X)));
    Y        = find(X > Xm);
   I4m(iv,ii)   = length(Y);
   I4a(iv,ii) = mean(X(Y));
end
I4t(ii)=length(POX6p);

end

XR1m = nanmean(R1m,1);
XR2m = nanmean(R2m,1);
XI3m = nanmean(I3m,1)./I3t;
XI4m = nanmean(I4m,1)./I4t;

XR1a= nanmean(R1a,1);
XR2a = nanmean(R2a,1);
XI3a = nanmean(I3a,1);
XI4a = nanmean(I4a,1);


[p,h,stat] = signrank(XI3m(SG),XI4m(SG))
[p,h,stat] = signrank(XI3m(AG),XI4m(AG))

[p,h,stat] = signrank(XI3a(SG),XI4a(SG))
[p,h,stat] = signrank(XI3a(AG),XI4a(AG))

writetable(table(XI3m(SG),XI4m(SG) ,XI3m(AG) ,XI4m(AG),'VariableNames',{'I3S','I4S', 'I3W','I4W'}), 'Bursts_Beta_desync_C3C4_Ratio_1.5s.txt');
writetable(table(XI3a(SG),XI4a(SG) ,XI3a(AG) ,XI4a(AG),'VariableNames',{'I3S','I4S', 'I3W','I4W'}), 'Bursts_Beta_desync_C3C4_Amp_1.5s.txt');


figure;
subplot(4,2,1); plot(POX5.time , POX5p(2,:)'); ylim([0 3])
subplot(4,2,2); plot(POX5.time , POX6p(4,:)'); ylim([0 3])
subplot(4,2,3); plot(POX5.time , POX5p(6,:)'); ylim([0 3])
subplot(4,2,4); plot(POX5.time , POX6p(7,:)'); ylim([0 3])
subplot(4,2,5); plot(POX5.time , POX5p(22,:)'); ylim([0 3])
subplot(4,2,6); plot(POX5.time , POX6p(40,:)'); ylim([0 3])
subplot(4,2,7); plot(POX5.time , POX5p(26,:)'); ylim([0 3])
subplot(4,2,8); plot(POX5.time , POX6p(17,:)'); ylim([0 3])


%% Correlation with Performance Behavior

%Performace Gain
PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);
PG_S = PG.Var1(SG); PG_W=PG.Var1(AG);

for i = 1:11
[r(i), p(i)] =corr(POS6_SB.powspctrm(:,i) , PG_S, 'type','Kendall');
% [r2(i), p2(i)] =corr(PRE6_SB.powspctrm(:,1) - PRE5_SB.powspctrm(:,1), C3_NE, 'type','Kendall')
end

for i = 1:11
[rr(i), pr(i)] =corr(PRE6_SB.powspctrm(:,i)  -  PRE5_SB.powspctrm(:,i), PG_S, 'type','kendall')

end
