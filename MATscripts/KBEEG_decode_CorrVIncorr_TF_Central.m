%% Start toolboxes

%fieldtrip
clear;
addpath(genpath('/home/b1044271/Toolboxes/obob_ownft'));
obob_init_ft;

%EEGLab
addpath(genpath('/home/b1044271/Toolboxes/eeglab14_1_1b'));
eeglab; close gcf;

%MVPA
addpath(genpath('/home/b1044271/Toolboxes/MVPA-Light'));
startup_MVPA_Light;
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

trials = NaN(33,2000);

%% classifier parameters
% Configuration struct for time classification with cross-validation. We
% perform 5-fold cross-validation with 2 repetitions. As classifier, we
% use LDA with its default settings.
cfg_LDA                 = [];
cfg_LDA.cv              = 'kfold';
cfg_LDA.k               = 5;
cfg_LDA.repeat          = 2;
cfg_LDA.classifier      = 'lda';
cfg_LDA.metric          = 'auc';
cfg_LDA.prerprocess     = 'zscore';
cfg_LDA.sample_dimension    = 1;
cfg_LDA.feature_dimension   = 2;
cfg_LDA.dimension_names = {'samples','channels','frequencies','time points'};

% cfg_LDA3=cfg_LDA;
% cfg_LDA3.hyperparameter          = [];
% cfg_LDA3.hyperparameter.lambda   = 'auto';
% cfg_sl = [];
% cfg_sl.average     = 1;
% cfg_sl.metric      = 'accuracy';
% cfg_sl.size        = 1; 

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
        EEG1 = pop_cleanline(EEG1, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG1.nbchan                            = EEG1.nbchan+1;
    EEG1.data(end+1,:)                     = zeros(1, EEG1.pnts);
    EEG1.chanlocs(1,EEG1.nbchan).labels    = 'CZ';
    EEG1                                   = pop_chanedit(EEG1, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG1                                   = eeg_checkset( EEG1);
    EEG1                                   = pop_reref( EEG1, []);  % normal eeglab reref    
    
        EEG2                                   = pop_eegfiltnew(EEG2, 0.1, []);
            EEG2 = pop_cleanline(EEG2, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

    EEG2.nbchan                            = EEG2.nbchan+1;
    EEG2.data(end+1,:)                     = zeros(1, EEG2.pnts);
    EEG2.chanlocs(1,EEG2.nbchan).labels    = 'CZ';
    EEG2                                   = pop_chanedit(EEG2, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG2                                   = eeg_checkset( EEG2);
    EEG2                                   = pop_reref( EEG2, []);  % normal eeglab reref    
    
%     EEG3.nbchan                          = EEG3.nbchan+1;
%     EEG3.data(end+1,:)                   = zeros(1, EEG3.pnts);
%     EEG3.chanlocs(1,EEG3.nbchan).labels  = 'CZ';
%     EEG3                                 = pop_chanedit(EEG3, 'lookup','C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
%     EEG3                                 = eeg_checkset( EEG3);
%     EEG3                                 = pop_reref( EEG3, []);  % normal eeglab reref    
%     EEG3                                 = pop_eegfiltnew(EEG3, 1, 60);
%     
%     EEG4.nbchan                          = EEG4.nbchan+1;
%     EEG4.data(end+1,:)                   = zeros(1, EEG4.pnts);
%     EEG4.chanlocs(1,EEG4.nbchan).labels  = 'CZ';
%     EEG4                                 = pop_chanedit(EEG4, 'lookup','C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
%     EEG4                                 = eeg_checkset( EEG4);
%     EEG4                                 = pop_reref( EEG4, []);  % normal eeglab reref    
%     EEG4                                 = pop_eegfiltnew(EEG4, 1, 60);
    
    EEG5                                 = pop_eegfiltnew(EEG5, 0.1, []);
    EEG5 = pop_cleanline(EEG5, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

    EEG5.nbchan                          = EEG5.nbchan+1;
    EEG5.data(end+1,:)                   = zeros(1, EEG5.pnts);
    EEG5.chanlocs(1,EEG5.nbchan).labels  = 'CZ';
    EEG5                                 = pop_chanedit(EEG5, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG5                                 = eeg_checkset( EEG5);
    EEG5                                 = pop_reref( EEG5, []);  % normal eeglab reref    
    
        EEG6                                 = pop_eegfiltnew(EEG6, 0.1, []);
    EEG6 = pop_cleanline(EEG6, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

    EEG6.nbchan                          = EEG6.nbchan+1;
    EEG6.data(end+1,:)                   = zeros(1, EEG6.pnts);
    EEG6.chanlocs(1,EEG6.nbchan).labels  = 'CZ';
    EEG6                                 = pop_chanedit(EEG6, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG6                                 = eeg_checkset( EEG6);
    EEG6                                 = pop_reref( EEG6, []);  % normal eeglab reref    
    
%     EEG7.nbchan                          = EEG7.nbchan+1;
%     EEG7.data(end+1,:)                   = zeros(1, EEG7.pnts);
%     EEG7.chanlocs(1,EEG7.nbchan).labels  = 'CZ';
%     EEG7                                 = pop_chanedit(EEG7, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
%     EEG7                                 = eeg_checkset( EEG7);
%     EEG7                                 = pop_reref( EEG7, []);  % normal eeglab reref    
%     EEG7                                 = pop_eegfiltnew(EEG7, 1, 60);

    % central
    EEG1 = pop_select(EEG1, 'channel', 4:5); EEG2 = pop_select(EEG2, 'channel', 4:5); 
%     EEG3 = pop_select(EEG3, 'channel', 1:10); EEG4 = pop_select(EEG4, 'channel', 1:10); 
    EEG5 = pop_select(EEG5, 'channel', 4:5); EEG6 = pop_select(EEG6, 'channel', 4:5); 
    
    %% Epoching
    % REG1    
%     R1_C = pop_epoch (EEG1,{'255'}, [-5 5]);
%     R1_I = pop_epoch (EEG1,{'254'}, [-5 5]);
   
    % REG2
%     R2_C = pop_epoch (EEG2,{'255'}, [-5 5]);
%     R2_I   = pop_epoch (EEG2,{'254'}, [-5 5]);

  % INV3
    I3_C = pop_epoch (EEG5,{'255'}, [-5 5]);
    I3_I = pop_epoch (EEG5,{'254'}, [-5 5]);

    % INV4
    I4_C = pop_epoch (EEG6,{'255'}, [-5 5]);
    I4_I = pop_epoch (EEG6,{'254'}, [-5 5]);
    
    trials(i,1) = length(I3_C.epoch); trials(i,2) = length(I3_I.epoch); trials(i,3) = length(I4_C.epoch); trials(i,4) = length(I4_I.epoch);
    
   %% epoching
    I3C2 =eeglab2fieldtrip(I3_C,'preprocessing','none');
    I3I2   =eeglab2fieldtrip(I3_I,'preprocessing','none');

    I4C2 =eeglab2fieldtrip(I4_C,'preprocessing','none');
    I4I2   =eeglab2fieldtrip(I4_I,'preprocessing','none');
    
     cfg=[];
  cfg.method      ='mtmconvol'; 
  cfg.taper       ='dpss';
  cfg.output      = 'pow';
  cfg.foi         = 00.5:0.5:45;
  cfg.t_ftimwin   = 5./cfg.foi;  % 5 cycles per time window
  cfg.toi         = -5:0.01:5;
  cfg.tapsmofrq   = 0.4*cfg.foi;  % 
  cfg.keeptrials  = 'yes'; % 

    
  TF_I3C{i} = ft_freqanalysis(cfg,I3C2);
  TF_I3I{i}   = ft_freqanalysis(cfg,I3I2);

  TF_I4C{i}=ft_freqanalysis(cfg,I4C2);
  TF_I4I{i}  =ft_freqanalysis(cfg,I4I2);
  
 %% Decoding
 %%%%%%%%%%
 X1 = cat(1, TF_I3C{i}.powspctrm, TF_I3I{i}.powspctrm);
 label_C = cat(1,ones(size(TF_I3C{i}.powspctrm,1),1), 2*ones(size(TF_I3I{i}.powspctrm,1),1));

 X2 = cat(1, TF_I4C{i}.powspctrm, TF_I4I{i}.powspctrm);
 label_I = cat(1,ones(size(TF_I4C{i}.powspctrm,1),1), 2*ones(size(TF_I4I{i}.powspctrm,1),1));

 %Run classification across time
[auc_I3{i}, result_I3{i}]     = mv_classify(cfg_LDA, X1, label_C);
% [aucRI_2{i}, resultRI_2{i}]     = mv_classify(cfg_LDA, X2, label_RI_2);
[auc_I4{i}, result_I4{i}]   = mv_classify(cfg_LDA, X2, label_I);

end


% MEAN OVER SUBJECTS
for i = 1:33    
    decode_I3(i,1:90,1:1001) = result_I3{i}.perf;
    decode_I4(i,1:90,1:1001)= result_I4{i}.perf;  
end

%% PERMUTATIONS

load('/home/b1044271/EEG_KB/tFV_all.mat')

eFV_all.powspctrm = permute(decode_I3,[1,4,2,3]);
eFV_all.time = TF_I3C{1,1}.time;
eFV_all.dimord = 'subj_chan_freq_time';
eFV_all.freq = TF_I3C{1,1}.freq;

I3_SG = eFV_all; I3_SG.powspctrm =  I3_SG.powspctrm(SG,:,:,:);
I3_WG = eFV_all; I3_WG.powspctrm =  I3_WG.powspctrm(AG,:,:,:);

eFV_all2 = eFV_all;
eFV_all2.powspctrm = permute(decode_I4,[1,4,2,3]);

I4_SG = eFV_all2; I4_SG.powspctrm =  I4_SG.powspctrm(SG,:,:,:);
I4_WG = eFV_all2; I4_WG.powspctrm =  I4_WG.powspctrm(AG,:,:,:);

cfg=[];
cfg.statistic        ='ft_statfun_depsamplesT'; %between groups independant sampels statistics
cfg.method           ='montecarlo'; % Monte Carlo method for calculating the signif prob an estimate of the p-value under the permutation distribution.
%  cfg.channel          = {'F3'};%frontal2
% cfg.channel          = {'C3','C4'};%frontal2
cfg.latency          = [-0.3 0.3];
 cfg.frequency        = [0.5 30];
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
  cfg.avgoverchan      = 'yes';
%cfg.minnbchan        = 2;
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

% SLEEP
[stat_f]   = ft_freqstatistics(cfg, I4_SG , I3_SG );

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

%Wake
[stat_f2]   = ft_freqstatistics(cfg, I4_WG , I3_WG );

%% COHEN D
cfg = [];
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.latency     = [-0.3 0.3];
cfg.avgoverchan = 'yes';
cfg.alpha       = 0.05;
cfg.numrandomization = 'all';
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.spmversion  = 'spm12';
cfg.uvar             = 1;
cfg.ivar             = 2;
cfg.parameter   = 'powspctrm';

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

inference1  = ft_freqstatistics(cfg,  I4_SG , I3_SG  );

%GRAND AVERAGE (optional)
cfg = [];
cfg.latency     = [-0.3 0.3];
grandavgFIC_sel = ft_selectdata(cfg,  I4_SG);
grandavgFIC_sel.powspctrm = squeeze(grandavgFIC_sel.powspctrm);


grandavgFC_sel  = ft_selectdata(cfg, I3_SG);
grandavgFC_sel.powspctrm = squeeze(grandavgFC_sel.powspctrm);

% COHEN'S D
x1 = nan(16,1);
x2 = nan(16,1);

for i=1:16

  %construct a 3-dimensional Boolean array to select the data from this participant
  sel3d = false(size(grandavgFIC_sel.powspctrm));
  sel3d(i,:,:) = inference1.negclusterslabelmat==1;

  %select the FIC data in the cluster for this participant, represent it as a vector
  tmp = grandavgFIC_sel.powspctrm(sel3d(:));
  
  %compute the average over the cluster
  x1(i) = mean(tmp);

  %select the FC data in the cluster for this participant, represent it as a vector
  tmp = grandavgFC_sel.powspctrm(sel3d(:));
  
  %compute the average over the cluster
  x2(i) = mean(tmp);
end

n1 = length(x1);
n2 = length(x2);

cd_TFKC= mean(x1-x2) ./ std(x1-x2)

