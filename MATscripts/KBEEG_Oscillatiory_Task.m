% KB EEG task 
%aperiodic activity KBEEG
%IRASA

% % start matlabpool with max local (physical) workers
% %-----------------------------------
% mypool = gcp;
% if mypool.Connected~=1
%     numCores = feature('numcores');
%     parpool(numCores);
% end

%% Start toolboxes
rng(20)
%fieldtrip
clear;
addpath(genpath('/home/b1044271/Toolboxes/obob_ownft'));
obob_init_ft;

%EEGLab
addpath(genpath('/home/b1044271/Toolboxes/eeglab14_1_1b'));
eeglab; close gcf;


%% load files

path_R='/home/b1044271/EEG_KB/AG_REG_Raw/';
dataRow_r   =struct2cell(dir(fullfile(path_R , '*.cnt'))');
files_R =dataRow_r(1,:);
%savepath= 'C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\AG_REG_Raw\preprocessed_normal_reref\';

path_I='/home/b1044271/EEG_KB/AG_INV_Raw/';
dataRow_i   =struct2cell(dir(fullfile(path_I , '*.cnt'))');
files_I =dataRow_i(1,:);

%% Groups 
AG = [1:10, 22, 23,24,25,31,32,33];
SG = [11:21, 26,27,28,29,30];


for i = 1:33

    if i <10%read data (replace with arnod delorms new function)
        EEG1=pop_loadcnt([path_R sprintf('VP0%dREG1.cnt', i)]);
        EEG2=pop_loadcnt([path_R sprintf('VP0%dREG2.cnt', i)]);
%         EEG3=pop_loadcnt([path_I sprintf('VP0%dINV1.cnt', i)]);
%         EEG4=pop_loadcnt([path_I sprintf('VP0%dINV2.cnt', i)]);
        EEG5=pop_loadcnt([path_I sprintf('VP0%dINV3.cnt', i)]);
        EEG6=pop_loadcnt([path_I sprintf('VP0%dINV4.cnt', i)]);
        EEG7=pop_loadcnt([path_I sprintf('VP0%dINV5.cnt', i)]);
    else
        EEG1=pop_loadcnt([path_R sprintf('VP%dREG1.cnt', i)]);
        EEG2=pop_loadcnt([path_R sprintf('VP%dREG2.cnt', i)]);
%         EEG3=pop_loadcnt([path_I sprintf('VP%dINV1.cnt', i)]);
%         EEG4=pop_loadcnt([path_I sprintf('VP%dINV2.cnt', i)]);
        EEG5=pop_loadcnt([path_I sprintf('VP%dINV3.cnt', i)]);
        EEG6=pop_loadcnt([path_I sprintf('VP%dINV4.cnt', i)]);
        EEG7=pop_loadcnt([path_I sprintf('VP%dINV5.cnt', i)]);
    end

    %% Prepro
    % delete mastoids, eye and muscles
    EEG1 = pop_select(EEG1, 'channel', 1:10); EEG2 = pop_select(EEG2, 'channel', 1:10); 
%     EEG3 = pop_select(EEG3, 'channel', 1:10); EEG4 = pop_select(EEG4, 'channel', 1:10); 
    EEG5 = pop_select(EEG5, 'channel', 1:10); EEG6 = pop_select(EEG6, 'channel', 1:10); 
    EEG7 = pop_select(EEG7, 'channel', 1:10); 
    
    %add the REF channel back and perform the average reference  
    EEG1                                   = pop_eegfiltnew(EEG1, 0.1, []); %pop_eegfiltnew(EEG1, 13, 30);
%     EEG1 = pop_cleanline(EEG1, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG1.nbchan                            = EEG1.nbchan+1;
    EEG1.data(end+1,:)                     = zeros(1, EEG1.pnts);
    EEG1.chanlocs(1,EEG1.nbchan).labels    = 'CZ';
    EEG1                                   = pop_chanedit(EEG1, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG1                                   = eeg_checkset( EEG1);
    EEG1                                   = pop_reref( EEG1, []);  % normal eeglab reref    
    
    EEG2                                   = pop_eegfiltnew(EEG2, 0.1, []);
%     EEG2 = pop_cleanline(EEG2, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG2.nbchan                            = EEG2.nbchan+1;
    EEG2.data(end+1,:)                     = zeros(1, EEG2.pnts);
    EEG2.chanlocs(1,EEG2.nbchan).labels    = 'CZ';
    EEG2                                   = pop_chanedit(EEG2, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG2                                   = eeg_checkset( EEG2);
    EEG2                                   = pop_reref( EEG2, []);  % normal eeglab reref    
 

    EEG5                                 = pop_eegfiltnew(EEG5, 0.1, []);
%     EEG5 = pop_cleanline(EEG5, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG5.nbchan                          = EEG5.nbchan+1;
    EEG5.data(end+1,:)                   = zeros(1, EEG5.pnts);
    EEG5.chanlocs(1,EEG5.nbchan).labels  = 'CZ';
    EEG5                                 = pop_chanedit(EEG5, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG5                                 = eeg_checkset( EEG5);
    EEG5                                 = pop_reref( EEG5, []);  % normal eeglab reref    

    EEG6                                 = pop_eegfiltnew(EEG6, 0.1, []);
%     EEG6 = pop_cleanline(EEG6, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG6.nbchan                          = EEG6.nbchan+1;
    EEG6.data(end+1,:)                   = zeros(1, EEG6.pnts);
    EEG6.chanlocs(1,EEG6.nbchan).labels  = 'CZ';
    EEG6                                 = pop_chanedit(EEG6, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG6                                 = eeg_checkset( EEG6);
    EEG6                                 = pop_reref( EEG6, []);  % normal eeglab reref    

%     EEG7                                 = pop_eegfiltnew(EEG7, 0.1, 60);
%     EEG7 = pop_cleanline(EEG7, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
%     EEG7.nbchan                          = EEG7.nbchan+1;
%     EEG7.data(end+1,:)                   = zeros(1, EEG7.pnts);
%     EEG7.chanlocs(1,EEG7.nbchan).labels  = 'CZ';
%     EEG7                                 = pop_chanedit(EEG7, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
%     EEG7                                 = eeg_checkset( EEG7);
%     EEG7                                 = pop_reref( EEG7, []);  % normal eeglab reref    
 
    %% epoching/segemnting
    EEG1 = eeg_regepochs (EEG1, 'recurrence', 5, 'limits', [0 10], 'eventtype', '1');
    EEG2 = eeg_regepochs (EEG2, 'recurrence', 5, 'limits', [0 10], 'eventtype', '1');
%     EEG3 = eeg_regepochs (EEG3, 'recurrence', 5, 'limits', [0 10], 'eventtype', '1');
%     EEG4 = eeg_regepochs (EEG4, 'recurrence', 5, 'limits', [0 10], 'eventtype', '1');
    EEG5 = eeg_regepochs (EEG5, 'recurrence', 5, 'limits', [0 10], 'eventtype', '1');
    EEG6 = eeg_regepochs (EEG6, 'recurrence', 5, 'limits', [0 10], 'eventtype', '1');
%     EEG7 = eeg_regepochs (EEG7, 'recurrence', 5, 'limits', [0 10], 'eventtype', '1');
    
    
    % run bad epoch detection and reject bad epoches
    [EEG1, rmepochs] = pop_autorej( EEG1,'nogui','on');epoch_props = epoch_properties(EEG1,1:EEG1.nbchan);
    BadEpochs   = min_z(epoch_props);EEG1 = pop_rejepoch(EEG1, BadEpochs,0);
    
    [EEG2, rmepochs] = pop_autorej( EEG2,'nogui','on');epoch_props = epoch_properties(EEG2,1:EEG2.nbchan);
    BadEpochs   = min_z(epoch_props);EEG2 = pop_rejepoch(EEG2, BadEpochs,0);

%     [EEG3, rmepochs] = pop_autorej( EEG3,'nogui','on');epoch_props = epoch_properties(EEG3,1:EEG3.nbchan);
%     BadEpochs   = min_z(epoch_props);EEG3 = pop_rejepoch(EEG3, BadEpochs,0);
%     
%     [EEG4, rmepochs] = pop_autorej( EEG4,'nogui','on');epoch_props = epoch_properties(EEG4,1:EEG4.nbchan);
%     BadEpochs   = min_z(epoch_props);EEG4 = pop_rejepoch(EEG4, BadEpochs,0);
    
    [EEG5, rmepochs] = pop_autorej( EEG5,'nogui','on');epoch_props = epoch_properties(EEG5,1:EEG5.nbchan);
    BadEpochs   = min_z(epoch_props);EEG5 = pop_rejepoch(EEG5, BadEpochs,0);
    
    [EEG6, rmepochs] = pop_autorej( EEG6,'nogui','on');epoch_props = epoch_properties(EEG6,1:EEG6.nbchan);
    BadEpochs   = min_z(epoch_props);EEG6 = pop_rejepoch(EEG6, BadEpochs,0);
    
%     [EEG7, rmepochs] = pop_autorej(EEG7,'nogui','on');epoch_props = epoch_properties(EEG7,1:EEG7.nbchan);
%     BadEpochs   = min_z(epoch_props);EEG7 = pop_rejepoch(EEG7, BadEpochs,0);
%     
    
    %% Equal number of epochs
%     minI = min([length(length(EEG5.epoch), length(EEG6.epoch)]); trials(i,1)=minI;
    minR = min([length(EEG1.epoch), length(EEG2.epoch), length(EEG5.epoch), length(EEG6.epoch)]); trials(i,2)=minR;

           
    EEG1e = pop_select(EEG1, 'trial', datasample(1:length(EEG1.epoch),minR,'Replace',false));
    EEG2e = pop_select(EEG2, 'trial', datasample(1:length(EEG2.epoch),minR,'Replace',false));
        
%     EEG3e = pop_select(EEG3, 'trial', datasample(1:length(EEG3.epoch),minI,'Replace',false));
%     EEG4e = pop_select(EEG4, 'trial', datasample(1:length(EEG4.epoch),minI,'Replace',false));

    EEG5e = pop_select(EEG5, 'trial', datasample(1:length(EEG5.epoch),minR,'Replace',false));
    EEG6e = pop_select(EEG6, 'trial', datasample(1:length(EEG6.epoch),minR,'Replace',false));
    
    
    %% IRASA
    EEG1f=eeglab2fieldtrip(EEG1e,'preprocessing','none');
    EEG2f=eeglab2fieldtrip(EEG2e,'preprocessing','none');
%     EEG3f=eeglab2fieldtrip(EEG3e,'preprocessing','none');
%     EEG4f=eeglab2fieldtrip(EEG4e,'preprocessing','none');
    EEG5f=eeglab2fieldtrip(EEG5e,'preprocessing','none');
    EEG6f=eeglab2fieldtrip(EEG6e,'preprocessing','none');
%     EEG7f=eeglab2fieldtrip(EEG7e,'preprocessing','none');
    
    cfg               = [];
    cfg.foi           = 1:0.5:90;
    cfg.pad           = 'nextpow2';
    cfg.method        = 'irasa';
    cfg.output        = 'fractal';
    cfg.keeptrials    = 'yes';
      
    EEG1a{i}=ft_freqanalysis(cfg,EEG1f);
    EEG2a{i}=ft_freqanalysis(cfg,EEG2f);
%     EEG3a{i}=ft_freqanalysis(cfg,EEG3f);
%     EEG4a{i}=ft_freqanalysis(cfg,EEG4f);
    EEG5a{i}=ft_freqanalysis(cfg,EEG5f);
    EEG6a{i}=ft_freqanalysis(cfg,EEG6f);
%     EEG7a{i}=ft_freqanalysis(cfg,EEG7f);
    
    cfg.output        = 'original';
    EEG1b{i}=ft_freqanalysis(cfg,EEG1f);
    EEG2b{i}=ft_freqanalysis(cfg,EEG2f);
%     EEG3b{i}=ft_freqanalysis(cfg,EEG3f);
%     EEG4b{i}=ft_freqanalysis(cfg,EEG4f);
    EEG5b{i}=ft_freqanalysis(cfg,EEG5f);
    EEG6b{i}=ft_freqanalysis(cfg,EEG6f);
%     EEG7b{i}=ft_freqanalysis(cfg,EEG7f);
    
    
    cfg                   = [];
    cfg.parameter         = 'powspctrm';
    cfg.operation         = 'x1-x2';
      
    EEG1o{i}=ft_math(cfg,EEG1b{i}, EEG1a{i});
    EEG2o{i}=ft_math(cfg,EEG2b{i}, EEG2a{i});
%     EEG3o{i}=ft_math(cfg,EEG3b{i}, EEG3a{i});
%     EEG4o{i}=ft_math(cfg,EEG4b{i}, EEG4a{i});
    EEG5o{i}=ft_math(cfg,EEG5b{i}, EEG5a{i});
    EEG6o{i}=ft_math(cfg,EEG6b{i}, EEG6a{i});
%     EEG7o{i}=ft_math(cfg,EEG7b{i}, EEG7a{i});
        
end    

% Descpritive
cfg_g =[];
for i = 1:33

EEG1x{i}=ft_freqdescriptives(cfg_g, EEG1o{i});
EEG2x{i}=ft_freqdescriptives(cfg_g, EEG2o{i});
EEG5x{i}=ft_freqdescriptives(cfg_g, EEG5o{i});
EEG6x{i}=ft_freqdescriptives(cfg_g, EEG6o{i});

EEG1fx{i}=ft_freqdescriptives(cfg_g, EEG1a{i});
EEG2fx{i}=ft_freqdescriptives(cfg_g, EEG2a{i});
EEG5fx{i}=ft_freqdescriptives(cfg_g, EEG5a{i});
EEG6fx{i}=ft_freqdescriptives(cfg_g, EEG6a{i});

end
%% GRAND AVERAGE of all subj N =33
cfg_g=[];
cfg_g.keepindividual = 'yes';

EEG1fr=ft_freqgrandaverage(cfg_g, EEG1fx{:}); 
EEG2fr=ft_freqgrandaverage(cfg_g, EEG2fx{:});
% EEG3fr=ft_freqgrandaverage(cfg_g, EEG3a{:});
% EEG4fr=ft_freqgrandaverage(cfg_g, EEG4a{:});
EEG5fr=ft_freqgrandaverage(cfg_g, EEG5fx{:});
EEG6fr=ft_freqgrandaverage(cfg_g, EEG6fx{:});
% EEG7fr=ft_freqgrandaverage(cfg_g, EEG7a{:});

%ALL OSCILLATORY
EEG1os=ft_freqgrandaverage(cfg_g, EEG1x{:});
EEG2os=ft_freqgrandaverage(cfg_g, EEG2x{:});
% EEG3os=ft_freqgrandaverage(cfg_g, EEG3o{:});
% EEG4os=ft_freqgrandaverage(cfg_g, EEG4o{:});
EEG5os=ft_freqgrandaverage(cfg_g, EEG5x{:});
EEG6os=ft_freqgrandaverage(cfg_g, EEG6x{:});
% EEG7os=ft_freqgrandaverage(cfg_g, EEG7o{:});

%ALL ORIG
% EEG1og=ft_freqgrandaverage(cfg_g, EEG1b{:});
% EEG2og=ft_freqgrandaverage(cfg_g, EEG2b{:});
% EEG3og=ft_freqgrandaverage(cfg_g, EEG3b{:});
% EEG4og=ft_freqgrandaverage(cfg_g, EEG4b{:});
% EEG5og=ft_freqgrandaverage(cfg_g, EEG5b{:});
% EEG6og=ft_freqgrandaverage(cfg_g, EEG6b{:});
% EEG7og=ft_freqgrandaverage(cfg_g, EEG7b{:});

%% SLOPE OF APERIODIC

cfg                            =[];
cfg.keepindividual     = 'yes';
cfg.frequency           = [1 45];
% cfg.channel           = {'C3','C4'};
% cfg.avgoverchan      = 'yes';

EEG1_frac = ft_selectdata(cfg, EEG1fr);
EEG2_frac = ft_selectdata(cfg, EEG2fr);
% EEG3_frac = ft_selectdata(cfg, EEG3fr);
% EEG4_frac = ft_selectdata(cfg, EEG4fr);
EEG5_frac = ft_selectdata(cfg, EEG5fr);
EEG6_frac = ft_selectdata(cfg, EEG6fr);
% EEG7_frac = ft_selectdata(cfg, EEG7fr);

for i = 1:33
    
pR1(i,:) = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG1_frac.powspctrm(i,:,:))), 1);
pR2(i,:) = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG2_frac.powspctrm(i,:,:))), 1); 
% pI1(i,:) = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG3_frac.powspctrm(i,:,:))), 1);
% pI2(i,:) = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG4_frac.powspctrm(i,:,:))), 1);  
pI3(i,:) = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG5_frac.powspctrm(i,:,:))), 1);
pI4(i,:) = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG6_frac.powspctrm(i,:,:))), 1);  
end

% ALL
sR1=pR1(:,1);
sR2=pR2(:,1);
% sI1=pI1(:,1);
% sI2=pI2(:,1);
sI3=pI3(:,1);
sI4=pI4(:,1);
% sI5=pI5(:,1);
writetable(table(sR1,sR2,sI3,sI4,'VariableNames',{'R1', 'R2','I3','I4'}), 'RvI_Slope145_C_X.txt')



for i = 1:33
     for ii = 1:11
pR1a = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG1_frac.powspctrm(i,ii,:))), 1); pR1(i,ii)=pR1a(:,1);
pR2a = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG2_frac.powspctrm(i,ii,:))), 1); pR2(i,ii)=pR1a(:,1);

pI3a = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG5_frac.powspctrm(i,ii,:))), 1); pI3(i,ii)=pI3a(:,1);
pI4a = polyfit(log(EEG1_frac.freq)',  squeeze(log(EEG6_frac.powspctrm(i,ii,:))), 1);  pI4(i,ii)=pI4a(:,1);
     end
end

pR1s = pR1(SG,:);
pR2s = pR2(SG,:);
pI3s = pI3(SG,:);
pI4s = pI4(SG,:);

%Performace Gain
PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);
PG_S = PG.Var1(SG); PG_W=PG.Var1(AG);

for i = 1:11
[r(i), p(i)] =corr(pI4s(:,i)- pI3s(:,i), PG_S, 'type','pearson')
end

% Topoplots

load('/home/b1044271/Toolboxes/LAYOUT/ChanInfoFT') % for my setup
idx = [ 1 2 3 4 11 5 6 7 8 9 10];
r1n = r(idx); p1n = p(idx);

% Create some dummy data
Data.label = lay.label; % Channel names
Data.dimord = 'chan_freq'; 
Data.freq = 1:10; %arbitrary frequency vector
Data.zvalues = repelem(r1n',1,length(Data.freq)); 

% Do topoplot 
cfg = [];
cfg.layout = lay; %layout struct
cfg.parameter = 'zvalues';
cfg.colormap=colormap(brewermap(256, '*RdYlBu'));   
cfg.style = 'straight';
cfg.gridscale =1000;
cfg.markersize =15;
cfg.highlightsize=60;
cfg.highlightcolor=[0 0 0];
cfg.highlightsymbol='.';
cfg.zlim = [-0.2 0.7];
cfg.highlight = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';

cfg.highlightchannel = find(p1n*11 < 0.05);
figure; ft_topoplotER(cfg,Data); title('Spectral slope (3-55Hz)'); 
set(gca,'FontSize',18,'TickLength',[.01 .01],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1.5); 
h = colorbar();
w = h.LineWidth;
h.LineWidth = 2;
ylabel(h, 't','FontSize', 22);
set(gca,'FontSize',18,'TickLength',[.01 .01],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1.5); 

X_pre = cat(1,sR1_sl, sI3_sl, sR1_w, sI3_w);
X_post = cat(1,sR2_sl, sI4_sl, sR2_w, sI4_w);
writetable(table(X_pre, X_post,'VariableNames',{'pre', 'post'}), 'RvI_Slope355.txt')

%% DELTA
cfg                =[];
cfg.keepindividual = 'yes';
cfg.frequency      = [0.1 4];
cfg.avgoverfreq    = 'yes';


EEG1_delta = ft_selectdata(cfg, EEG1os);
EEG2_delta = ft_selectdata(cfg, EEG2os);
EEG3_delta = ft_selectdata(cfg, EEG3os);
EEG4_delta = ft_selectdata(cfg, EEG4os);
EEG5_delta = ft_selectdata(cfg, EEG5os);
EEG6_delta = ft_selectdata(cfg, EEG6os);
EEG7_delta = ft_selectdata(cfg, EEG7os);

%% OSCILLATIONS - THETA (4-8 Hz)

cfg                        = [];
cfg.keepindividual = 'yes';
cfg.frequency       = [4 8];
cfg.avgoverfreq   = 'yes';
% cfg.channel           = {'Pz','P4'};
cfg.avgoverchan   = 'yes';

EEG1_theta = ft_selectdata(cfg, EEG1os);
EEG2_theta = ft_selectdata(cfg, EEG2os);
% EEG3_theta = ft_selectdata(cfg, EEG3os);
% EEG4_theta = ft_selectdata(cfg, EEG4os);
EEG5_theta = ft_selectdata(cfg, EEG5os);
EEG6_theta = ft_selectdata(cfg, EEG6os);
EEG7_theta = ft_selectdata(cfg, EEG7os);

PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);
PG_S = PG.Var1(SG); PG_W=PG.Var1(AG);

for i = 1:11
[r(i), p(i)] =corr(EEG6_theta.powspctrm(SG,i)-EEG5_theta.powspctrm(SG,i), PG_S, 'type','Kendall')
end

%% OSCILLATIONS - Alpha (8-12 Hz)

cfg                =[];
cfg.keepindividual = 'yes';
cfg.frequency      = [8 12];
cfg.avgoverfreq    = 'yes';
% cfg.channel           = {'Pz','P4'};
cfg.avgoverchan   = 'yes';

EEG1_alpha = ft_selectdata(cfg, EEG1os);
EEG2_alpha = ft_selectdata(cfg, EEG2os);
% EEG3_alpha = ft_selectdata(cfg, EEG3os);
% EEG4_alpha = ft_selectdata(cfg, EEG4os);
EEG5_alpha = ft_selectdata(cfg, EEG5os);
EEG6_alpha = ft_selectdata(cfg, EEG6os);
% EEG7_alpha = ft_selectdata(cfg, EEG7os);

%% OSCILLATIONS - BETA (13-30 Hz)

cfg                =[];
cfg.keepindividual = 'yes';
cfg.frequency      = [13 30];
cfg.avgoverfreq    = 'yes';
cfg.channel            = {'C3','C4'};
  cfg.avgoverchan    = 'yes';

EEG1_beta = ft_selectdata(cfg, EEG1os);
EEG2_beta = ft_selectdata(cfg, EEG2os);
% EEG3_beta = ft_selectdata(cfg, EEG3os);
% EEG4_beta = ft_selectdata(cfg, EEG4os);
EEG5_beta = ft_selectdata(cfg, EEG5os);
EEG6_beta = ft_selectdata(cfg, EEG6os);
% EEG7_beta = ft_selectdata(cfg, EEG7os);

for i = 1:11
[r(i), p(i)] =corr(EEG6_beta.powspctrm(SG,i)-EEG5_beta.powspctrm(SG,i), PG_S, 'type','Kendall')
end

%% KB PLots

% APERIOIDC
f=figure;subplot(3,2,1);
SG_F = mean([sI1_sl, sI2_sl,sI3_sl,sI4_sl,sI5_sl, sR1_sl,sR2_sl]);
WG_F = mean([sI1_w, sI2_w,sI3_w,sI4_w,sI5_w, sR1_w,sR2_w]);

eSG_F = std([sI1_sl, sI2_sl,sI3_sl,sI4_sl,sI5_sl, sR1_sl, sR2_sl ])./sqrt(16);
eWG_F = std([sI1_w, sI2_w,sI3_w,sI4_w,sI5_w, sR1_w, sR2_w])./sqrt(17);
colors1 = [1 0 0; 0 0 1; 0 1 0];
bl(:,1)=boundedline([1,2,3,4,5,7,8], SG_F, eSG_F ,[1,2,3,4,5,7,8], WG_F, eWG_F,'cmap', colors1, 'alpha','transparency', 0.2);

% Delta
subplot(3,2,2);
SG_F = mean([EEG3_delta.powspctrm(SG), EEG4_delta.powspctrm(SG),EEG5_delta.powspctrm(SG),EEG6_delta.powspctrm(SG),EEG7_delta.powspctrm(SG), EEG1_delta.powspctrm(SG),EEG2_delta.powspctrm(SG)]);
WG_F = mean([EEG3_delta.powspctrm(AG), EEG4_delta.powspctrm(AG),EEG5_delta.powspctrm(AG),EEG6_delta.powspctrm(AG),EEG7_delta.powspctrm(AG), EEG1_delta.powspctrm(AG),EEG2_delta.powspctrm(AG)]);

eSG_F = std([EEG3_delta.powspctrm(SG), EEG4_delta.powspctrm(SG),EEG5_delta.powspctrm(SG),EEG6_delta.powspctrm(SG),EEG7_delta.powspctrm(SG), EEG1_delta.powspctrm(SG),EEG2_delta.powspctrm(SG)])./sqrt(16);
eWG_F = std([EEG3_delta.powspctrm(AG), EEG4_delta.powspctrm(AG),EEG5_delta.powspctrm(AG),EEG6_delta.powspctrm(AG),EEG7_delta.powspctrm(AG), EEG1_delta.powspctrm(AG),EEG2_delta.powspctrm(AG)])./sqrt(17);
colors1 = [1 0 0; 0 0 1; 0 1 0];
bl(:,2)=boundedline([1,2,3,4,5,7,8], SG_F, eSG_F ,[1,2,3,4,5,7,8], WG_F, eWG_F,'cmap', colors1, 'alpha','transparency', 0.2);

% Theta
subplot(3,2,3);
SG_F = mean([EEG3_theta.powspctrm(SG), EEG4_theta.powspctrm(SG),EEG5_theta.powspctrm(SG),EEG6_theta.powspctrm(SG),EEG7_theta.powspctrm(SG), EEG1_theta.powspctrm(SG),EEG2_theta.powspctrm(SG)]);
WG_F = mean([EEG3_theta.powspctrm(AG), EEG4_theta.powspctrm(AG),EEG5_theta.powspctrm(AG),EEG6_theta.powspctrm(AG),EEG7_theta.powspctrm(AG), EEG1_theta.powspctrm(AG),EEG2_theta.powspctrm(AG)]);

eSG_F = std([EEG3_theta.powspctrm(SG), EEG4_theta.powspctrm(SG),EEG5_theta.powspctrm(SG),EEG6_theta.powspctrm(SG),EEG7_theta.powspctrm(SG), EEG1_theta.powspctrm(SG),EEG2_theta.powspctrm(SG)])./sqrt(16);
eWG_F = std([EEG3_theta.powspctrm(AG), EEG4_theta.powspctrm(AG),EEG5_theta.powspctrm(AG),EEG6_theta.powspctrm(AG),EEG7_theta.powspctrm(AG), EEG1_theta.powspctrm(AG),EEG2_theta.powspctrm(AG)])./sqrt(17);
colors1 = [1 0 0; 0 0 1; 0 1 0];
bl(:,3)=boundedline([1,2,3,4,5,7,8], SG_F, eSG_F ,[1,2,3,4,5,7,8], WG_F, eWG_F,'cmap', colors1, 'alpha','transparency', 0.2);

% Alpha
subplot(3,2,4);
SG_F = mean([EEG3_alpha.powspctrm(SG), EEG4_alpha.powspctrm(SG),EEG5_alpha.powspctrm(SG),EEG6_alpha.powspctrm(SG),EEG7_alpha.powspctrm(SG), EEG1_alpha.powspctrm(SG),EEG2_alpha.powspctrm(SG)]);
WG_F = mean([EEG3_alpha.powspctrm(AG), EEG4_alpha.powspctrm(AG),EEG5_alpha.powspctrm(AG),EEG6_alpha.powspctrm(AG),EEG7_alpha.powspctrm(AG), EEG1_alpha.powspctrm(AG),EEG2_alpha.powspctrm(AG)]);

eSG_F = std([EEG3_alpha.powspctrm(SG), EEG4_alpha.powspctrm(SG),EEG5_alpha.powspctrm(SG),EEG6_alpha.powspctrm(SG),EEG7_alpha.powspctrm(SG), EEG1_alpha.powspctrm(SG),EEG2_alpha.powspctrm(SG)])./sqrt(16);
eWG_F = std([EEG3_alpha.powspctrm(AG), EEG4_alpha.powspctrm(AG),EEG5_alpha.powspctrm(AG),EEG6_alpha.powspctrm(AG),EEG7_alpha.powspctrm(AG), EEG1_alpha.powspctrm(AG),EEG2_alpha.powspctrm(AG)])./sqrt(17);
colors1 = [1 0 0; 0 0 1; 0 1 0];
bl(:,4)=boundedline([1,2,3,4,5,7,8], SG_F, eSG_F ,[1,2,3,4,5,7,8], WG_F, eWG_F,'cmap', colors1, 'alpha','transparency', 0.2);

% BETA
subplot(3,2,5);
SG_F = mean([EEG3_beta.powspctrm(SG), EEG4_beta.powspctrm(SG),EEG5_beta.powspctrm(SG),EEG6_beta.powspctrm(SG),EEG7_beta.powspctrm(SG), EEG1_beta.powspctrm(SG),EEG2_beta.powspctrm(SG)]);
WG_F = mean([EEG3_beta.powspctrm(AG), EEG4_beta.powspctrm(AG),EEG5_beta.powspctrm(AG),EEG6_beta.powspctrm(AG),EEG7_beta.powspctrm(AG), EEG1_beta.powspctrm(AG),EEG2_beta.powspctrm(AG)]);

eSG_F = std([EEG3_beta.powspctrm(SG), EEG4_beta.powspctrm(SG),EEG5_beta.powspctrm(SG),EEG6_beta.powspctrm(SG),EEG7_beta.powspctrm(SG), EEG1_beta.powspctrm(SG),EEG2_beta.powspctrm(SG)])./sqrt(16);
eWG_F = std([EEG3_beta.powspctrm(AG), EEG4_beta.powspctrm(AG),EEG5_beta.powspctrm(AG),EEG6_beta.powspctrm(AG),EEG7_beta.powspctrm(AG), EEG1_beta.powspctrm(AG),EEG2_beta.powspctrm(AG)])./sqrt(17);
colors1 = [1 0 0; 0 0 1; 0 1 0];
bl(:,5)=boundedline([1,2,3,4,5,7,8], SG_F, eSG_F ,[1,2,3,4,5,7,8], WG_F, eWG_F,'cmap', colors1, 'alpha','transparency', 0.2);

set(bl(:),'LineStyle','-','LineWidth',4) %set line width
lh = legend(bl(:)); title('Beta');
legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';
ax = findobj(f,'Type','Axes');
namess = {'Aperiodic','Delta','Theta','Alpha', 'Beta'};
Xs = {'I1','I2','I3','I4', 'I5'};
for i=1:length(ax)
%      hXLabel(i) = xlabel('\fontsize{26} Session');
%      ylim(ax(i),[-7 1])
    xlim(ax(i),[1 5])
    set(ax(i),'XTick', 1:1:5,'XTickLabel', Xs);
    title(ax(i),namess{6-i})
    set(ax(i),'FontName','sans serif','FontSize',14);
%   set([ hXLabel(i) hYLabel(i) ], 'FontName','sans serif','FontSize',24);
 set(ax(i),'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
end

han=axes(f,'visible','off');
han.XLabel.Visible='on';
xlabel(han,sprintf('\n Sessions'), 'FontName','sans serif','FontSize',24);
han.YLabel.Visible='on';
ylabel(han,{'Power (uv2)', ' '}, 'FontName','sans serif','FontSize',24);

%% PERCENTAGE CHANGE KB PLOTS

figure();
SG_F = mean([((sI2_sl-sI1_sl)./sI1_sl)*100,((sI3_sl-sI1_sl)./sI1_sl)*100,((sI4_sl-sI1_sl)./sI1_sl)*100,((sI5_sl-sI1_sl)./sI1_sl)*100]);
WG_F = mean([((sI2_w-sI1_w)./sI1_w)*100,((sI3_w-sI1_w)./sI1_w)*100,((sI4_w-sI1_w)./sI1_w)*100,((sI5_w-sI1_w)./sI1_w)*100]);

eSG_F = std([((sI2_sl-sI1_sl)./sI1_sl)*100,((sI3_sl-sI1_sl)./sI1_sl)*100,((sI4_sl-sI1_sl)./sI1_sl)*100,((sI5_sl-sI1_sl)./sI1_sl)*100])./sqrt(16);
eWG_F = std([((sI2_w-sI1_w)./sI1_w)*100,((sI3_w-sI1_w)./sI1_w)*100,((sI4_w-sI1_w)./sI1_w)*100,((sI5_w-sI1_w)./sI1_w)*100])./sqrt(17);

colors1 = [0 0 1; 1 0 0];
bl = boundedline([1,2,3,4,5], [0 SG_F], [0 eSG_F] ,[1,2,3,4,5], [0 WG_F], [0 eWG_F],'cmap', colors1, 'alpha','transparency', 0.2);
line([3.5 3.5], [-40 50], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');
lh = legend(bl); set(bl,'LineStyle','-','LineWidth',4) %set line width
legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'southwest';
hXLabel = xlabel('\fontsize{26} Session');
hYLabel = ylabel('\fontsize{26} Slope change (%)');
set(gca,'FontName','sans serif','FontSize',14);
set([ hXLabel, hYLabel lh],'FontName','sans serif','FontSize',24);
ylim([-10 20])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
Xs = {'I1','I2','I3','I4', 'I5'};
set(gca,'XTick', 1:1:5,'XTickLabel', Xs, 'YTick', [-20, 0,20]);
title(gca,'Aperiodic (1/f)');

f=figure;
axesHandles(1) = subplot(2,2,1);
% Delta
SG_F = mean([ ((EEG4_delta.powspctrm(SG)-EEG3_delta.powspctrm(SG))./EEG3_delta.powspctrm(SG))*100,((EEG5_delta.powspctrm(SG)-EEG3_delta.powspctrm(SG))./EEG3_delta.powspctrm(SG))*100,...
    ((EEG6_delta.powspctrm(SG)-EEG3_delta.powspctrm(SG))./EEG3_delta.powspctrm(SG))*100, ((EEG7_delta.powspctrm(SG)-EEG3_delta.powspctrm(SG))./EEG3_delta.powspctrm(SG))*100]);
    
WG_F = mean([ ((EEG4_delta.powspctrm(AG)-EEG3_delta.powspctrm(AG))./EEG3_delta.powspctrm(AG))*100,((EEG5_delta.powspctrm(AG)-EEG3_delta.powspctrm(AG))./EEG3_delta.powspctrm(AG))*100,...
    ((EEG6_delta.powspctrm(AG)-EEG3_delta.powspctrm(AG))./EEG3_delta.powspctrm(AG))*100, ((EEG7_delta.powspctrm(AG)-EEG3_delta.powspctrm(AG))./EEG3_delta.powspctrm(AG))*100]);


eSG_F = std([ ((EEG4_delta.powspctrm(SG)-EEG3_delta.powspctrm(SG))./EEG3_delta.powspctrm(SG))*100,((EEG5_delta.powspctrm(SG)-EEG3_delta.powspctrm(SG))./EEG3_delta.powspctrm(SG))*100,...
    ((EEG6_delta.powspctrm(SG)-EEG3_delta.powspctrm(SG))./EEG3_delta.powspctrm(SG))*100, ((EEG7_delta.powspctrm(SG)-EEG3_delta.powspctrm(SG))./EEG3_delta.powspctrm(SG))*100])./sqrt(16);
    
eWG_F = std([ ((EEG4_delta.powspctrm(AG)-EEG3_delta.powspctrm(AG))./EEG3_delta.powspctrm(AG))*100,((EEG5_delta.powspctrm(AG)-EEG3_delta.powspctrm(AG))./EEG3_delta.powspctrm(AG))*100,...
    ((EEG6_delta.powspctrm(AG)-EEG3_delta.powspctrm(AG))./EEG3_delta.powspctrm(AG))*100, ((EEG7_delta.powspctrm(AG)-EEG3_delta.powspctrm(AG))./EEG3_delta.powspctrm(AG))*100])./sqrt(16);

bl(:,1) = boundedline([1,2,3,4,5], [0 SG_F], [0 eSG_F] ,[1,2,3,4,5], [0 WG_F], [0 eWG_F],'cmap', colors1, 'alpha','transparency', 0.2);
line([3.5 3.5], [-40 50], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');

% Theta
axesHandles(2) = subplot(2,2,2);
SG_F = mean([ ((EEG4_theta.powspctrm(SG)-EEG3_theta.powspctrm(SG))./EEG3_theta.powspctrm(SG))*100,((EEG5_theta.powspctrm(SG)-EEG3_theta.powspctrm(SG))./EEG3_theta.powspctrm(SG))*100,...
    ((EEG6_theta.powspctrm(SG)-EEG3_theta.powspctrm(SG))./EEG3_theta.powspctrm(SG))*100, ((EEG7_theta.powspctrm(SG)-EEG3_theta.powspctrm(SG))./EEG3_theta.powspctrm(SG))*100]);
    
WG_F = mean([ ((EEG4_theta.powspctrm(AG)-EEG3_theta.powspctrm(AG))./EEG3_theta.powspctrm(AG))*100,((EEG5_theta.powspctrm(AG)-EEG3_theta.powspctrm(AG))./EEG3_theta.powspctrm(AG))*100,...
    ((EEG6_theta.powspctrm(AG)-EEG3_theta.powspctrm(AG))./EEG3_theta.powspctrm(AG))*100, ((EEG7_theta.powspctrm(AG)-EEG3_theta.powspctrm(AG))./EEG3_theta.powspctrm(AG))*100]);


eSG_F = std([ ((EEG4_theta.powspctrm(SG)-EEG3_theta.powspctrm(SG))./EEG3_theta.powspctrm(SG))*100,((EEG5_theta.powspctrm(SG)-EEG3_theta.powspctrm(SG))./EEG3_theta.powspctrm(SG))*100,...
    ((EEG6_theta.powspctrm(SG)-EEG3_theta.powspctrm(SG))./EEG3_theta.powspctrm(SG))*100, ((EEG7_theta.powspctrm(SG)-EEG3_theta.powspctrm(SG))./EEG3_theta.powspctrm(SG))*100])./sqrt(16);
    
eWG_F = std([ ((EEG4_theta.powspctrm(AG)-EEG3_theta.powspctrm(AG))./EEG3_theta.powspctrm(AG))*100,((EEG5_theta.powspctrm(AG)-EEG3_theta.powspctrm(AG))./EEG3_theta.powspctrm(AG))*100,...
    ((EEG6_theta.powspctrm(AG)-EEG3_theta.powspctrm(AG))./EEG3_theta.powspctrm(AG))*100, ((EEG7_theta.powspctrm(AG)-EEG3_theta.powspctrm(AG))./EEG3_theta.powspctrm(AG))*100])./sqrt(16);

bl(:,2) = boundedline([1,2,3,4,5], [0 SG_F], [0 eSG_F] ,[1,2,3,4,5], [0 WG_F], [0 eWG_F],'cmap', colors1, 'alpha','transparency', 0.2);
line([3.5 3.5], [-40 50], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');
% ylabel('\fontsize{26} Power change (%)');

% ALPHA
axesHandles(3) = subplot(2,2,3);
SG_F = mean([ ((EEG4_alpha.powspctrm(SG)-EEG3_alpha.powspctrm(SG))./EEG3_alpha.powspctrm(SG))*100,((EEG5_alpha.powspctrm(SG)-EEG3_alpha.powspctrm(SG))./EEG3_alpha.powspctrm(SG))*100,...
    ((EEG6_alpha.powspctrm(SG)-EEG3_alpha.powspctrm(SG))./EEG3_alpha.powspctrm(SG))*100, ((EEG7_alpha.powspctrm(SG)-EEG3_alpha.powspctrm(SG))./EEG3_alpha.powspctrm(SG))*100]);
    
WG_F = mean([ ((EEG4_alpha.powspctrm(AG)-EEG3_alpha.powspctrm(AG))./EEG3_alpha.powspctrm(AG))*100,((EEG5_alpha.powspctrm(AG)-EEG3_alpha.powspctrm(AG))./EEG3_alpha.powspctrm(AG))*100,...
    ((EEG6_alpha.powspctrm(AG)-EEG3_alpha.powspctrm(AG))./EEG3_alpha.powspctrm(AG))*100, ((EEG7_alpha.powspctrm(AG)-EEG3_alpha.powspctrm(AG))./EEG3_alpha.powspctrm(AG))*100]);


eSG_F = std([ ((EEG4_alpha.powspctrm(SG)-EEG3_alpha.powspctrm(SG))./EEG3_alpha.powspctrm(SG))*100,((EEG5_alpha.powspctrm(SG)-EEG3_alpha.powspctrm(SG))./EEG3_alpha.powspctrm(SG))*100,...
    ((EEG6_alpha.powspctrm(SG)-EEG3_alpha.powspctrm(SG))./EEG3_alpha.powspctrm(SG))*100, ((EEG7_alpha.powspctrm(SG)-EEG3_alpha.powspctrm(SG))./EEG3_alpha.powspctrm(SG))*100])./sqrt(16);
    
eWG_F = std([ ((EEG4_alpha.powspctrm(AG)-EEG3_alpha.powspctrm(AG))./EEG3_alpha.powspctrm(AG))*100,((EEG5_alpha.powspctrm(AG)-EEG3_alpha.powspctrm(AG))./EEG3_alpha.powspctrm(AG))*100,...
    ((EEG6_alpha.powspctrm(AG)-EEG3_alpha.powspctrm(AG))./EEG3_alpha.powspctrm(AG))*100, ((EEG7_alpha.powspctrm(AG)-EEG3_alpha.powspctrm(AG))./EEG3_alpha.powspctrm(AG))*100])./sqrt(16);

bl(:,3) = boundedline([1,2,3,4,5], [0 SG_F], [0 eSG_F] ,[1,2,3,4,5], [0 WG_F], [0 eWG_F],'cmap', colors1, 'alpha','transparency', 0.2);
line([3.5 3.5], [-40 50], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');

% BETA
axesHandles(4) = subplot(2,2,4);
SG_F = mean([ ((EEG4_beta.powspctrm(SG)-EEG3_beta.powspctrm(SG))./EEG3_beta.powspctrm(SG))*100,((EEG5_beta.powspctrm(SG)-EEG3_beta.powspctrm(SG))./EEG3_beta.powspctrm(SG))*100,...
    ((EEG6_beta.powspctrm(SG)-EEG3_beta.powspctrm(SG))./EEG3_beta.powspctrm(SG))*100, ((EEG7_beta.powspctrm(SG)-EEG3_beta.powspctrm(SG))./EEG3_beta.powspctrm(SG))*100]);
    
WG_F = mean([ ((EEG4_beta.powspctrm(AG)-EEG3_beta.powspctrm(AG))./EEG3_beta.powspctrm(AG))*100,((EEG5_beta.powspctrm(AG)-EEG3_beta.powspctrm(AG))./EEG3_beta.powspctrm(AG))*100,...
    ((EEG6_beta.powspctrm(AG)-EEG3_beta.powspctrm(AG))./EEG3_beta.powspctrm(AG))*100, ((EEG7_beta.powspctrm(AG)-EEG3_beta.powspctrm(AG))./EEG3_beta.powspctrm(AG))*100]);

eSG_F = std([ ((EEG4_beta.powspctrm(SG)-EEG3_beta.powspctrm(SG))./EEG3_beta.powspctrm(SG))*100,((EEG5_beta.powspctrm(SG)-EEG3_beta.powspctrm(SG))./EEG3_beta.powspctrm(SG))*100,...
    ((EEG6_beta.powspctrm(SG)-EEG3_beta.powspctrm(SG))./EEG3_beta.powspctrm(SG))*100, ((EEG7_beta.powspctrm(SG)-EEG3_beta.powspctrm(SG))./EEG3_beta.powspctrm(SG))*100])./sqrt(16);
    
eWG_F = std([ ((EEG4_beta.powspctrm(AG)-EEG3_beta.powspctrm(AG))./EEG3_beta.powspctrm(AG))*100,((EEG5_beta.powspctrm(AG)-EEG3_beta.powspctrm(AG))./EEG3_beta.powspctrm(AG))*100,...
    ((EEG6_beta.powspctrm(AG)-EEG3_beta.powspctrm(AG))./EEG3_beta.powspctrm(AG))*100, ((EEG7_beta.powspctrm(AG)-EEG3_beta.powspctrm(AG))./EEG3_beta.powspctrm(AG))*100])./sqrt(16);

bl(:,4) = boundedline([1,2,3,4,5], [0 SG_F], [0 eSG_F] ,[1,2,3,4,5], [0 WG_F], [0 eWG_F],'cmap', colors1, 'alpha','transparency', 0.2);
line([3.5 3.5], [-40 50], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');

set(bl(:),'LineStyle','-','LineWidth',4) %set line width
lh = legend(bl(:)); title('Beta');
legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';
ax = findobj(f,'Type','Axes');
namess = {'Aperiodic','Delta','Theta','Alpha', 'Beta'};
Xs = {'I1','I2','I3','I4', 'I5'};
for i=1:length(ax)
%      hXLabel(i) = xlabel('\fontsize{26} Session');
    ylim(ax(i),[-30 50])
    xlim(ax(i),[1 5])
    set(ax(i),'XTick', 1:1:5,'XTickLabel', Xs, 'YTick', [-30, 0,40]);
    title(ax(i),namess{6-i})
    set(ax(i),'FontName','sans serif','FontSize',14);
%   set([ hXLabel(i) hYLabel(i) ], 'FontName','sans serif','FontSize',24);
 set(ax(i),'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
end

han=axes(f,'visible','off');
han.XLabel.Visible='on';
xlabel(han,sprintf('\n Sessions'), 'FontName','sans serif','FontSize',24);
han.YLabel.Visible='on';
ylabel(han,{'Power change (%)', ' '}, 'FontName','sans serif','FontSize',24);

%% COrrelate performance gains with change in beta

beta_diff = EEG6_beta.powspctrm(SG) -  EEG5_beta.powspctrm(SG);
P_gain = table2array(readtable('/home/b1044271/EEG_KB/Results/Behavioral/PgainS.txt', 'Headerlines', 0));

[p,h]=corr(pp, beta_diff, 'Type','pearson');

%% REG - Power change

% APERIOIDC
figure;
SG_F = mean([sR1_sl, sR2_sl]);
WG_F = mean([sR1_w, sR2_w]);

colors1 = [0 0 1; 1 0 0; 0 1 0];
bl=boxplot([sR1_sl; sR2_sl;sR1_w; sR2_w], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ]);
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [sR1_sl; sR2_sl;sR1_w; sR2_w],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
lh = legend(bl(:));
legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';
hXLabel = xlabel('\fontsize{26} Session');
hYLabel = ylabel('\fontsize{26} Slope change (%)');
set(gca,'FontName','sans serif','FontSize',14);
set([ hXLabel, hYLabel lh],'FontName','sans serif','FontSize',24);
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

f=figure;
% Delta
subplot(2,2,1);
bl=boxplot([EEG1_delta.powspctrm(SG);EEG2_delta.powspctrm(SG);EEG1_delta.powspctrm(AG); EEG2_delta.powspctrm(AG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ]);
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [EEG1_delta.powspctrm(SG); EEG2_delta.powspctrm(SG);EEG1_delta.powspctrm(AG); EEG2_delta.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
set(gca,'FontName','sans serif','FontSize',14);
lh = legend(bl(:));
legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';
set( lh,'FontName','sans serif','FontSize',24);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

% Theta
subplot(2,2,2);
boxplot([EEG1_theta.powspctrm(SG);EEG2_theta.powspctrm(SG);EEG1_theta.powspctrm(AG); EEG2_theta.powspctrm(AG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ])
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [EEG1_theta.powspctrm(SG); EEG2_theta.powspctrm(SG);EEG1_theta.powspctrm(AG); EEG2_theta.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 


% Alpha
subplot(2,2,3);
boxplot([EEG1_alpha.powspctrm(SG);EEG2_alpha.powspctrm(SG);EEG1_alpha.powspctrm(AG); EEG2_alpha.powspctrm(AG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ])
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [EEG1_alpha.powspctrm(SG); EEG2_alpha.powspctrm(SG);EEG1_alpha.powspctrm(AG); EEG2_alpha.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 


% BETA
subplot(2,2,4);
boxplot([EEG1_beta.powspctrm(SG);EEG2_beta.powspctrm(SG);EEG1_beta.powspctrm(AG); EEG2_beta.powspctrm(AG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ])
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [EEG1_beta.powspctrm(SG); EEG2_beta.powspctrm(SG);EEG1_beta.powspctrm(AG); EEG2_beta.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
set(gca,'FontName','sans serif','FontSize',14);
set([ hXLabel, hYLabel lh],'FontName','sans serif','FontSize',24);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

ax = findobj(f,'Type','Axes');
namess = {'Delta','Theta','Alpha', 'Beta'};
for i=1:length(ax)
%      hXLabel(i) = xlabel('\fontsize{26} Session');
    title(ax(i),namess{5-i})
    set(ax(i),'FontName','sans serif','FontSize',14);
    
%   set([ hXLabel(i) hYLabel(i) ], 'FontName','sans serif','FontSize',24);
end

han=axes(f,'visible','off');
han.XLabel.Visible='on';
xlabel(han,sprintf('\n Sessions'), 'FontName','sans serif','FontSize',24);
han.YLabel.Visible='on';
ylabel(han,{'Power change (%)', ' '}, 'FontName','sans serif','FontSize',24);


%% REG - Power change

% APERIOIDC
figure;
SG_F = mean([s_sl, sR2_sl]);
WG_F = mean([sR1_w, sR2_w]);

colors1 = [0 0 1; 1 0 0; 0 1 0];
bl=boxplot([sR1_sl; sR2_sl;sR1_w; sR2_w], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ]);
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [sR1_sl; sR2_sl;sR1_w; sR2_w],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
lh = legend(bl(:));
legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';
hXLabel = xlabel('\fontsize{26} Session');
hYLabel = ylabel('\fontsize{26} Slope change (%)');
set(gca,'FontName','sans serif','FontSize',14);
set([ hXLabel, hYLabel lh],'FontName','sans serif','FontSize',24);
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

f=figure;
% Delta
subplot(2,2,1);
bl=boxplot([EEG1_delta.powspctrm(SG);EEG2_delta.powspctrm(SG);EEG1_delta.powspctrm(AG); EEG2_delta.powspctrm(AG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ]);
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [EEG1_delta.powspctrm(SG); EEG2_delta.powspctrm(SG);EEG1_delta.powspctrm(AG); EEG2_delta.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
set(gca,'FontName','sans serif','FontSize',14);
lh = legend(bl(:));
legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';
set( lh,'FontName','sans serif','FontSize',24);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

% Theta
subplot(2,2,2);
boxplot([EEG1_theta.powspctrm(SG);EEG2_theta.powspctrm(SG);EEG1_theta.powspctrm(AG); EEG2_theta.powspctrm(AG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ])
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [EEG1_theta.powspctrm(SG); EEG2_theta.powspctrm(SG);EEG1_theta.powspctrm(AG); EEG2_theta.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 


% Alpha
subplot(2,2,3);
boxplot([EEG1_alpha.powspctrm(SG);EEG2_alpha.powspctrm(SG);EEG1_alpha.powspctrm(AG); EEG2_alpha.powspctrm(AG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'b' 'b' ])
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1)), [EEG1_alpha.powspctrm(SG); EEG2_alpha.powspctrm(SG);EEG1_alpha.powspctrm(AG); EEG2_alpha.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

% BETA
f=figure;subplot(1,2,1);
boxplot([EEG1_beta.powspctrm(SG);EEG2_beta.powspctrm(SG);EEG3_beta.powspctrm(SG);EEG6_beta.powspctrm(SG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)),'Colors', ['b' 'b' 'b' 'b' ])
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)), [EEG1_beta.powspctrm(SG);EEG2_beta.powspctrm(SG);EEG3_beta.powspctrm(SG);EEG6_beta.powspctrm(SG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'I1 Sleep', 'I4 Sleep'})
set(gca,'FontName','sans serif','FontSize',14);
set([ hXLabel, hYLabel lh],'FontName','sans serif','FontSize',24);
ylim([0 0.02])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 


subplot(1,2,2);boxplot([EEG1_beta.powspctrm(AG);EEG2_beta.powspctrm(AG);EEG3_beta.powspctrm(AG) ;EEG6_beta.powspctrm(AG) ], cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'r' 'r' ])
hold on;scatter(cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)),[EEG1_beta.powspctrm(AG);EEG2_beta.powspctrm(AG);EEG3_beta.powspctrm(AG) ;EEG6_beta.powspctrm(AG) ],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 Wake', 'R2 Wake', 'I1 Wake', 'I4 Wake'})
set(gca,'FontName','sans serif','FontSize',14);
set([ hXLabel, hYLabel lh],'FontName','sans serif','FontSize',24);
 ylim([0 .02])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
% 
% ax = findobj(f,'Type','Axes');
% namess = {'Delta','Theta','Alpha', 'Beta'};
% for i=1:length(ax)
% %      hXLabel(i) = xlabel('\fontsize{26} Session');
%     title(ax(i),namess{5-i})
%     set(ax(i),'FontName','sans serif','FontSize',14);
%     
% %   set([ hXLabel(i) hYLabel(i) ], 'FontName','sans serif','FontSize',24);
% end

han=axes(f,'visible','off');
han.XLabel.Visible='on';
xlabel(han,sprintf('\n Sessions'), 'FontName','sans serif','FontSize',26);
han.YLabel.Visible='on';
ylabel(han,{'Power uv2', ' ',' '}, 'FontName','sans serif','FontSize',26);


%% APERIOIDC
figure;
SG_F = mean([sR1_sl, sR2_sl, sI3_sl, sI4_sl]);
WG_F = mean([sR1_w, sR2_w , sI3_w, sI4_w]);

colors1 = [0 0 1; 1 0 0; 0 1 0];
bl=boxplot([sR1_sl; sR2_sl; sR1_w; sR2_w; sI3_sl; sI4_sl;  sI3_w; sI4_w], cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),...
4*ones(17,1), 5*ones(16,1),6*ones(16,1), 7*ones(17,1),8*ones(17,1)),'Colors', ['b' 'b' 'r' 'r'  'b' 'b' 'r' 'r']);

hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(17,1),4*ones(17,1), 5*ones(16,1),6*ones(16,1), 7*ones(17,1),8*ones(17,1)),...
[sR1_sl; sR2_sl; sR1_w; sR2_w; sI3_sl; sI4_sl; sI3_w; sI4_w],60, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'R2 sleep', 'R1 Wake', 'R2 Wake', 'I3 sleep','I4 sleep', 'I3 Wake', 'I4 Wake'})
lh = legend(bl(:));
legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})]; 
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';
hXLabel = xlabel('\fontsize{26} Session');
hYLabel = ylabel('\fontsize{26} Slope change (%)');
set(gca,'FontName','sans serif','FontSize',14);
set([ hXLabel, hYLabel lh],'FontName','sans serif','FontSize',24);
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 


writetable(table(sR1,sR2,sI3,sI4, 'VariableNames',{'R1','R2','I3','I4'}),'FracSlopes_R_I_3_55.txt')

%time resolved
%ALL OSCILLATORY
EEG1os=ft_freqgrandaverage(cfg_g, EEG1o{:});
EEG2os=ft_freqgrandaverage(cfg_g, EEG2o{:});
EEG3os=ft_freqgrandaverage(cfg_g, EEG3o{:});
EEG4os=ft_freqgrandaverage(cfg_g, EEG4o{:});
EEG5os=ft_freqgrandaverage(cfg_g, EEG5o{:});
EEG6os=ft_freqgrandaverage(cfg_g, EEG6o{:});
EEG7os=ft_freqgrandaverage(cfg_g, EEG7o{:});

%% R vs I 

% SLEEP
%%%%%%%%%%

f=figure;
% Delta
subplot(2,2,1);
bl=boxplot([EEG1_delta.powspctrm(SG);EEG5_delta.powspctrm(SG);EEG2_delta.powspctrm(SG); EEG6_delta.powspctrm(SG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)),'Colors', ['b' 'b' 'b' 'b' ]);
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)), [EEG1_delta.powspctrm(SG);EEG5_delta.powspctrm(SG);EEG2_delta.powspctrm(SG); EEG6_delta.powspctrm(SG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'I3 sleep','R2 sleep', 'I4 sleep'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

% Theta
subplot(2,2,2);
bl=boxplot([EEG1_theta.powspctrm(SG);EEG5_theta.powspctrm(SG);EEG2_theta.powspctrm(SG); EEG6_theta.powspctrm(SG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)),'Colors', ['b' 'b' 'b' 'b' ]);
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)), [EEG1_theta.powspctrm(SG);EEG5_theta.powspctrm(SG);EEG2_theta.powspctrm(SG); EEG6_theta.powspctrm(SG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'I3 sleep','R2 sleep', 'I4 sleep'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

% Alpha
subplot(2,2,3);
bl=boxplot([EEG1_alpha.powspctrm(SG);EEG5_alpha.powspctrm(SG);EEG2_alpha.powspctrm(SG); EEG6_alpha.powspctrm(SG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)),'Colors', ['b' 'b' 'b' 'b' ]);
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)), [EEG1_alpha.powspctrm(SG);EEG5_alpha.powspctrm(SG);EEG2_alpha.powspctrm(SG); EEG6_alpha.powspctrm(SG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'I3 sleep','R2 sleep', 'I4 sleep'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

%BETA
subplot(2,2,4);
bl=boxplot([EEG1_beta.powspctrm(SG);EEG5_beta.powspctrm(SG);EEG2_beta.powspctrm(SG); EEG6_beta.powspctrm(SG)], cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)),'Colors', ['b' 'b' 'b' 'b' ]);
hold on;scatter(cat(1,ones(16,1),2*ones(16,1), 3*ones(16,1),4*ones(16,1)), [EEG1_beta.powspctrm(SG);EEG5_beta.powspctrm(SG);EEG2_beta.powspctrm(SG); EEG6_beta.powspctrm(SG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 sleep', 'I3 sleep','R2 sleep', 'I4 sleep'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

% WAKE
%%%%%%%%%%%%5
f=figure;
% Delta
subplot(2,2,1);
bl=boxplot([EEG1_delta.powspctrm(AG);EEG5_delta.powspctrm(AG);EEG2_delta.powspctrm(AG); EEG6_delta.powspctrm(AG)], cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'r' 'r' ]);
hold on;scatter(cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)), [EEG1_delta.powspctrm(AG);EEG5_delta.powspctrm(AG);EEG2_delta.powspctrm(AG); EEG6_delta.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 wake', 'I3 wake','R2 wake', 'I4 wake'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

% Theta
subplot(2,2,2);
bl=boxplot([EEG1_theta.powspctrm(AG);EEG5_theta.powspctrm(AG);EEG2_theta.powspctrm(AG); EEG6_theta.powspctrm(AG)], cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'r' 'r' ]);
hold on;scatter(cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)), [EEG1_theta.powspctrm(AG);EEG5_theta.powspctrm(AG);EEG2_theta.powspctrm(AG); EEG6_theta.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 wake', 'I3 wake','R2 wake', 'I4 wake'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

% Alpha
subplot(2,2,3);
bl=boxplot([EEG1_alpha.powspctrm(AG);EEG5_alpha.powspctrm(AG);EEG2_alpha.powspctrm(AG); EEG6_alpha.powspctrm(AG)], cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'r' 'r' ]);
hold on;scatter(cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)), [EEG1_alpha.powspctrm(AG);EEG5_alpha.powspctrm(AG);EEG2_alpha.powspctrm(AG); EEG6_alpha.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 wake', 'I3 wake','R2 wake', 'I4 wake'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

%BETA
subplot(2,2,4);
bl=boxplot([EEG1_beta.powspctrm(AG);EEG5_beta.powspctrm(AG);EEG2_beta.powspctrm(AG); EEG6_beta.powspctrm(AG)], cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)),'Colors', ['r' 'r' 'r' 'r' ]);
hold on;scatter(cat(1,ones(17,1),2*ones(17,1), 3*ones(17,1),4*ones(17,1)), [EEG1_beta.powspctrm(AG);EEG5_beta.powspctrm(AG);EEG2_beta.powspctrm(AG); EEG6_beta.powspctrm(AG)],50, [0 0 0],'filled')
set(findobj(gca,'type','line'),'linew',4) 
xticklabels({'R1 wake', 'I3 wake','R2 wake', 'I4 wake'})
set(gca,'FontName','sans serif','FontSize',14);
% ylim([0 1.5])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

%BETA
X_pre = cat(1,EEG1_beta.powspctrm(SG), EEG5_beta.powspctrm(SG), EEG1_beta.powspctrm(AG) ,EEG5_beta.powspctrm(AG));
X_post = cat(1,EEG2_beta.powspctrm(SG), EEG6_beta.powspctrm(SG), EEG2_beta.powspctrm(AG) ,EEG6_beta.powspctrm(AG));
writetable(table(X_pre, X_post,'VariableNames',{'pre', 'post'}), 'RvI_beta_C_X2.txt')

%Delta
X_pre = cat(1,EEG1_delta.powspctrm(SG), EEG5_delta.powspctrm(SG), EEG1_delta.powspctrm(AG) ,EEG5_delta.powspctrm(AG));
X_post = cat(1,EEG2_delta.powspctrm(SG), EEG6_delta.powspctrm(SG), EEG2_delta.powspctrm(AG) ,EEG6_delta.powspctrm(AG));
writetable(table(X_pre, X_post,'VariableNames',{'pre', 'post'}), 'RvI_delta.txt')

%Theta
X_pre = cat(1,EEG1_theta.powspctrm(SG), EEG5_theta.powspctrm(SG), EEG1_theta.powspctrm(AG) ,EEG5_theta.powspctrm(AG));
X_post = cat(1,EEG2_theta.powspctrm(SG), EEG6_theta.powspctrm(SG), EEG2_theta.powspctrm(AG) ,EEG6_theta.powspctrm(AG));
writetable(table(X_pre, X_post,'VariableNames',{'pre', 'post'}), 'RvI_theta_all_X.txt')


%alpha
X_pre = cat(1,EEG1_alpha.powspctrm(SG), EEG5_alpha.powspctrm(SG), EEG1_alpha.powspctrm(AG) ,EEG5_alpha.powspctrm(AG));
X_post = cat(1,EEG2_alpha.powspctrm(SG), EEG6_alpha.powspctrm(SG), EEG2_alpha.powspctrm(AG) ,EEG6_alpha.powspctrm(AG));
writetable(table(X_pre, X_post,'VariableNames',{'pre', 'post'}), 'RvI_alpha_all_X.txt')