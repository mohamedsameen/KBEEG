% KB EEG task 
%preprocessing and decoding 
%prepro: filter, remove channels, rereferencing and restoring Cz, and epoching
%Decode: Remove fractal component and decode correct and incorrect trials
%(LDA using MVPA light) on the oscillatory component


%% Start toolboxes

%fieldtrip
clear;
addpath(genpath('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\obob_ownft'));
obob_init_ft;

%EEGLab
addpath(genpath('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\eeglab14_1_1b'));
eeglab; close gcf;

%MVPA
addpath(genpath('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\MVPA-Light'));
startup_MVPA_Light;

% Colormap
addpath(genpath('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\DrosteEffect-BrewerMap-5b84f95'));
%% load files
cd 'C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\AG_REG_Raw'
path_R='C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\AG_REG_Raw\';
dataRow_r   =struct2cell(dir(fullfile(path_R , '*.cnt'))');
files_R =dataRow_r(1,:);
savepath= 'C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\AG_REG_Raw\preprocessed_normal_reref\';

cd 'C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\AG_INV_Raw'
path_I='C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\AG_INV_Raw\';
dataRow_i   =struct2cell(dir(fullfile(path_I , '*.cnt'))');
files_I =dataRow_i(1,:);

%% Groups 
AG = [1:10, 22, 23,24,25,31,32,33];
SG = [11:21, 26,27,28,29,30];

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
cfg_LDA.feature_dimension   = 3;        % now the time points act as features
cfg_LDA.dimension_names     = {'samples' 'electrodes' 'time points'};


for i = 1:33

    if i <10%read data (replace with arnod delorms new function)
        EEG1=pop_loadcnt([path_R sprintf('VP0%dREG1.cnt', i)]);
        EEG2=pop_loadcnt([path_R sprintf('VP0%dREG2.cnt', i)]);
%         EEG3=pop_loadcnt([path_I sprintf('VP0%dINV1.cnt', i)]);
%         EEG4=pop_loadcnt([path_I sprintf('VP0%dINV2.cnt', i)]);
        EEG5=pop_loadcnt([path_I sprintf('VP0%dINV3.cnt', i)]);
        EEG6=pop_loadcnt([path_I sprintf('VP0%dINV4.cnt', i)]);
%         EEG7=pop_loadcnt([path_I sprintf('VP0%dINV5.cnt', i)]);
    else
        EEG1=pop_loadcnt([path_R sprintf('VP%dREG1.cnt', i)]);
        EEG2=pop_loadcnt([path_R sprintf('VP%dREG2.cnt', i)]);
% %         EEG3=pop_loadcnt([path_I sprintf('VP%dINV1.cnt', i)]);
% %         EEG4=pop_loadcnt([path_I sprintf('VP%dINV2.cnt', i)]);
        EEG5=pop_loadcnt([path_I sprintf('VP%dINV3.cnt', i)]);
        EEG6=pop_loadcnt([path_I sprintf('VP%dINV4.cnt', i)]);
%         EEG7=pop_loadcnt([path_I sprintf('VP%dINV5.cnt', i)]);
    end

    %% Prepro
    % delete mastoids, eye and muscles
    EEG1 = pop_select(EEG1, 'channel', 1:10); EEG2 = pop_select(EEG2, 'channel', 1:10); 
%     EEG3 = pop_select(EEG3, 'channel', 1:10); EEG4 = pop_select(EEG4, 'channel', 1:10); 
    EEG5 = pop_select(EEG5, 'channel', 1:10); EEG6 = pop_select(EEG6, 'channel', 1:10); 
%     EEG7 = pop_select(EEG7, 'channel', 1:10); 
    
    % downsample to 125 Hz for decoding
%     EEG1 = pop_resample(EEG1,125); EEG2 = pop_resample(EEG2,125); 
%     EEG3 = pop_resample(EEG3,125); EEG4 = pop_resample(EEG4,125); 
%     EEG5 = pop_resample(EEG5,125); EEG6 = pop_resample(EEG6,125); 
%     EEG7 = pop_resample(EEG7,125);
%     
    %add the REF channel back and perform the average reference 
    EEG1                                   = pop_eegfiltnew(EEG1, 1, 45); %pop_eegfiltnew(EEG1, 13, 30);
%     EEG1 = pop_resample(EEG1,125);
%     EEG1 = pop_cleanline(EEG1, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG1.nbchan                            = EEG1.nbchan+1;
    EEG1.data(end+1,:)                     = zeros(1, EEG1.pnts);
    EEG1.chanlocs(1,EEG1.nbchan).labels    = 'CZ';
    EEG1                                   = pop_chanedit(EEG1, 'lookup','C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG1                                   = eeg_checkset( EEG1);
    EEG1                                   = pop_reref( EEG1, []);  % normal eeglab reref    
   
    EEG2                                   = pop_eegfiltnew(EEG2,1, 45);
% EEG2 = pop_resample(EEG2,125);   
%     EEG2 = pop_cleanline(EEG2, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
EEG2.nbchan                            = EEG2.nbchan+1;
    EEG2.data(end+1,:)                     = zeros(1, EEG2.pnts);
    EEG2.chanlocs(1,EEG2.nbchan).labels    = 'CZ';
    EEG2                                   = pop_chanedit(EEG2, 'lookup','C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG2                                   = eeg_checkset( EEG2);
    EEG2                                   = pop_reref( EEG2, []);  % normal eeglab reref    
 
    EEG5                                 = pop_eegfiltnew(EEG5, 1, 45);
%         EEG5 = pop_cleanline(EEG5, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
% EEG5 = pop_resample(EEG5,125);   
EEG5.nbchan                          = EEG5.nbchan+1;
    EEG5.data(end+1,:)                   = zeros(1, EEG5.pnts);
    EEG5.chanlocs(1,EEG5.nbchan).labels  = 'CZ';
    EEG5                                 = pop_chanedit(EEG5, 'lookup','C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG5                                 = eeg_checkset( EEG5);
    EEG5                                 = pop_reref( EEG5, []);  % normal eeglab reref    
   
    EEG6                                 = pop_eegfiltnew(EEG6, 1, 45);
%     EEG6 = pop_resample(EEG6,125);
%     EEG6 = pop_cleanline(EEG6, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG6.nbchan                          = EEG6.nbchan+1;
    EEG6.data(end+1,:)                   = zeros(1, EEG6.pnts);
    EEG6.chanlocs(1,EEG6.nbchan).labels  = 'CZ';
    EEG6                                 = pop_chanedit(EEG6, 'lookup','C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEG6                                 = eeg_checkset( EEG6);
    EEG6                                 = pop_reref( EEG6, []);  % normal eeglab reref    
  
%     EEG7                                 = pop_eegfiltnew(EEG7, 0.1, 45);
%     EEG7 = pop_cleanline(EEG7, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
%     EEG7.nbchan                          = EEG7.nbchan+1;
%     EEG7.data(end+1,:)                   = zeros(1, EEG7.pnts);
%     EEG7.chanlocs(1,EEG7.nbchan).labels  = 'CZ';
%     EEG7                                 = pop_chanedit(EEG7, 'lookup','C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
%     EEG7                                 = eeg_checkset( EEG7);
%     EEG7                                 = pop_reref( EEG7, []);  % normal eeglab reref    
 
    %% Epoching
    % REG1    
    R1_letters = pop_epoch (EEG1,{'255','254'}, [-0.3 0.3]);

    % REG2
    R2_letters = pop_epoch (EEG2,{'255','254'}, [-0.3 0.3]);

    % INV3
    I3_letters = pop_epoch (EEG5,{'255','254'}, [-0.3 0.3]);


    % INV4
    I4_letters = pop_epoch (EEG6,{'255','254'}, [-0.3 0.3]);


    % INV5
%     I5_letters = pop_epoch (EEG7,{'255','254'}, [-0.3 0.3]);

    
    %% Equal epochs
    
%     min_pre = min([length(R1_letters.epoch) length(I3_letters.epoch)]); 
%     min_post = min([length(R2_letters.epoch) length(I4_letters.epoch)]);
%     minT(i,1) = min_pre;
%     minT(i,2) = min_post;
    
%     R1_new = pop_select(R1_letters, 'trial', datasample(1:length(R1_letters.epoch),min_pre,'Replace',false));
%     I3_new = pop_select(I3_letters, 'trial', datasample(1:length(I3_letters.epoch),min_pre,'Replace',false));
%     R2_new = pop_select(R2_letters, 'trial', datasample(1:length(R2_letters.epoch),min_post,'Replace',false));
%     I4_new = pop_select(I4_letters, 'trial', datasample(1:length(I4_letters.epoch),min_post,'Replace',false));
% 
% 
%     min_post2 = min([length(I5_letters.epoch) length(I4_letters.epoch)]);
%     I4_new2   = pop_select(I4_letters, 'trial', datasample(1:length(I4_letters.epoch),min_post2,'Replace',false));
%     I5_new    = pop_select(I5_letters, 'trial', datasample(1:length(I5_letters.epoch),min_post2,'Replace',false));

    %% convert to FT
    preR1=eeglab2fieldtrip(R1_letters,'preprocessing','none');
    preR2=eeglab2fieldtrip(R2_letters,'preprocessing','none');
    preI3=eeglab2fieldtrip(I3_letters,'preprocessing','none');
    preI4=eeglab2fieldtrip(I4_letters,'preprocessing','none');
%     preI5=eeglab2fieldtrip(I5_letters,'preprocessing','none');

    %% Time-lock
  cfg=[];
  cfg.keeptrials = 'yes';
  ERP_preR1{i}=ft_timelockanalysis(cfg,preR1);
  ERP_preR2{i}=ft_timelockanalysis(cfg,preR2);
  ERP_preI3{i}=ft_timelockanalysis(cfg,preI3); 
  ERP_preI4{i}=ft_timelockanalysis(cfg,preI4);
%   ERP_preI5{i}=ft_timelockanalysis(cfg,preI5);

%% Baseline 
cfg_b=[];
cfg_b.baseline=[-0.3 -0.1];
cfg_b.baselinetype='relchange';
  
ERP_preR12{i} = ft_timelockbaseline(cfg_b,ERP_preR1{i});
ERP_preR22{i} = ft_timelockbaseline(cfg_b,ERP_preR2{i});
ERP_preI32{i} = ft_timelockbaseline(cfg_b,ERP_preI3{i});
ERP_preI42{i} = ft_timelockbaseline(cfg_b,ERP_preI4{i});

  %% concatenate
  
  %before retention
  X1 = cat(1, ERP_preR12{i}.trial, ERP_preI32{i}.trial);
   label_RI_1 = cat(1,ones(size(ERP_preR12{i}.trial,1),1), 2*ones(size(ERP_preI32{i}.trial,1),1));

  %after retention
  X2 = cat(1, ERP_preR22{i}.trial, ERP_preI42{i}.trial);
  label_RI_2 = cat(1,ones(size(ERP_preR22{i}.trial,1),1), 2*ones(size(ERP_preI42{i}.trial,1),1));

 
  %% Decode across time
  %Run classification across time

[aucRI_1{i}, resultRI_1{i}]     = mv_classify_across_time(cfg_LDA, X1, label_RI_1);
[aucRI_2{i}, resultRI_2{i}]     = mv_classify_across_time(cfg_LDA, X2, label_RI_2);

[aucRI_1S{i}, resultRI_1S{i}]     = mv_classify(cfg_LDA, X1, label_RI_1);
[aucRI_2S{i}, resultRI_2S{i}]     = mv_classify(cfg_LDA, X2, label_RI_2);

end

for i = 1:33
    decode_RI1(1:300,i) = resultRI_1{i}.perf;
    decode_RI2(1:300,i) = resultRI_2{i}.perf;
%     decode_I45(1:300,i) = resultI45{i}.perf;
    
    decode_RI1S(1:11,i) = resultRI_1S{i}.perf;
    decode_RI2S(1:11,i) = resultRI_2S{i}.perf;
%     decode_I45S(1:11,i) = resultI45_S{i}.perf;

end

decodeRI1_SG = decode_RI1(:,SG);
decodeRI1_WG = decode_RI1(:,AG);
decodeRI2_SG = decode_RI2(:,SG);
decodeRI2_WG = decode_RI2(:,AG);
% decodeI45_SG = decode_I45(:,SG);
% decodeI45_WG = decode_I45(:,AG);

decodeRI1S_SG = decode_RI1S(:,SG);
decodeRI1S_WG = decode_RI1S(:,AG);
decodeRI2S_SG = decode_RI2S(:,SG);
decodeRI2S_WG = decode_RI2S(:,AG);
% decodeI45S_SG = decode_I45S(:,SG);
% decodeI45S_WG = decode_I45S(:,AG);

% writetable(table(decodeRI1_SG),'Decode_RI1_SG.txt', 'WriteVariableNames', 0)
% writetable(table(decodeRI2_SG),'Decode_RI2_SG.txt', 'WriteVariableNames', 0)
% writetable(table(decodeRI1_WG),'Decode_RI1_WG.txt', 'WriteVariableNames', 0)
% writetable(table(decodeRI2_WG),'Decode_RI1_WG.txt', 'WriteVariableNames', 0)
% save('decodeRI1_SG', 'decodeRI1_SG');
% save('decodeRI2_SG', 'decodeRI2_SG');
% save('decodeRI1_WG', 'decodeRI1_WG');
% save('decodeRI2_WG', 'decodeRI2_WG');

%% TOPOPLOTS
%sleep
temp          = resultRI_1S{1};
temp.perf     = mean(decodeRI2S_SG,2) - mean(decodeRI1S_SG,2);
temp.perf_std = std(decodeRI1S_SG,0,2);
mv_plot_result(temp)
%wake
temp2 = temp;
temp2.perf     = mean(decodeRI2S_WG,2) - mean(decodeRI1S_WG,2);
temp2.perf_std = std(decodeRI2S_WG,0,2);

%adjust channel order
P(1:4,1)  = temp.perf(1:4);P2(1:4,1)  = temp2.perf(1:4);
P(5,1)    = temp.perf(11);P2(5,1)    = temp2.perf(11);
P(6:11,1) = temp.perf(5:10); P2(6:11,1) = temp2.perf(5:10);

[V,I]=maxk(abs(P),2); %max 3 electrodes
[V2,I2]=maxk(abs(P2),2);



ChM = zeros(11,1);
ChM(I) = 1;
load('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\LAYOUT\ChanInfoFT.mat')
cfg_plot = [];
cfg_plot.outline = lay.outline;
cfg_plot.res=2000;
cfg_plot.mark_chans = ChM;
cfg_plot.cbtitle = 'auc (Post - Pre)';
cfg_plot.markersize =10;
cfg_plot.clim = [-0.05 0.05];
figure;
mv_plot_topography2(cfg_plot, P, lay.pos);
colormap(brewermap(256, '*RdYlBu'))

%%%%%%%%%%%%%%%%%%


cfg.highlightchannel = find(ismember(p2n,maxk(p2n,3)));
subplot(1,2,2); ft_topoplotER(cfg,Data2); title('After Sleep'); 
set(gca,'FontSize',28,'TickLength',[.02 .02],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5); 

cfg.zlim = [-0.1 0.1];
cfg.highlightchannel = find(ismember(abs(Data3.zvalues),maxk(abs(Data3.zvalues),3)));
figure; ft_topoplotER(cfg,Data3); title('Difference'); 
set(gca,'FontSize',28,'TickLength',[.02 .02],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5); 

%% Collapsed stats

% [p1,h1,stat1]=signrank(RI2a_WG, RI1a_WG)
% [p2,h2,stat2]=signrank(RI2a_SG, RI1a_SG)
% [p3,h3,stat3]=signrank(RI2aa_SG, RI1a_SG) %Only I4
% [p3a,h3a,stat3a]=signrank(RI2aa_WG, RI1a_WG) %only I4
% 
% [p4,h4,stat4]=signrank(RI2b_WG, RI1b_WG)
% [p5,h5,stat5]=signrank(RI2b_SG, RI1b_SG)
% 
% [p6,h6,stat6]=ranksum(RI2ab_SG, RI2ab_WG)

%% Within Subject Permutation
%%%%%%%%%%%%%%%%%
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
cfg.design= design;  % design matrix

% design the fieldtrip structures
cfg_g=[];
cfg_g.channel = 'C3';
cfg_g.keepindividual = 'yes';
pre_R2S=ft_timelockgrandaverage(cfg_g, ERP_preR2{SG});

Before_S = pre_R2S; Before_S.individual =  permute(decodeRI1_SG,[2,3,1]);
After_S  = pre_R2S; After_S.individual =  permute(decodeRI2_SG,[2,3,1]);


% run stats
[stat_S]   = ft_timelockstatistics(cfg, After_S, Before_S);

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

Before_W = pre_R2S; Before_W.individual =  permute(decodeRI1_WG,[2,3,1]);
After_W  = pre_R2S; After_W.individual =  permute(decodeRI2_WG,[2,3,1]);
% run stats
[stat_W]   = ft_timelockstatistics(cfg, After_W, Before_W);

%% BETWEEN SUBJECT DESIGN 
cfg=[];
cfg.statistic        ='ft_statfun_indepsamplesT'; %between groups independant sampels statistics
cfg.method           ='montecarlo'; % Monte Carlo method for calculating the signif prob an estimate of the p-value under the permutation distribution.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.avgoverchan      = 'yes';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
cfg.ivar             = 1;

design = [ones(1,16), ones(1,17)*2];
cfg.design= design;  % design matrix

% [stat_SW2]   = ft_timelockstatistics(cfg, After_S, After_W);
% [stat_SW]   = ft_timelockstatistics(cfg, Before_S, Before_W);

Diff_S = After_S; Diff_S.individual = After_S.individual - Before_S.individual;
Diff_W = After_W; Diff_W.individual = After_W.individual - Before_W.individual;

[stat_diff]   = ft_timelockstatistics(cfg, Diff_S, Diff_W);

%% Plotting time series data 
% %decode across time 
 f = figure;
% Sleep - Sleep
subplot(2,1,1);colors1 = [0 0.5 0 ;0.25 0.25 0.25 ];
bl=boundedline(ERP_preI3{i}.time,mean(decodeRI2_SG,2), std(decodeRI2_SG')./sqrt(16),ERP_preI3{i}.time,mean(decodeRI1_SG,2), std(decodeRI1_SG')./sqrt(16),'cmap', colors1, 'alpha','transparency', 0.2);
set(bl,'LineStyle','-','LineWidth',4) %set line width
line([0 0], [-250 250], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');
%line([-0.07 0.27], [0.7 0.7], 'color','k','LineStyle','-', 'LineWidth',3, 'DisplayName', ''); text(-0.04,0.72,'* p = 0.019','FontName','sans serif','FontSize',20);
%  line([-0.2 -0.13], [0.7 0.7], 'color','k','LineStyle','-', 'LineWidth',3, 'DisplayName', '');text(-0.2,0.72,'* p = 0.005','FontName','sans serif','FontSize',20);
 line([-0.3 -0.258], [0.65 0.65], 'color','k','LineStyle','-', 'LineWidth',4, 'DisplayName', '');text(-0.3,0.67,'* p = 0.03','FontName','sans serif','FontSize',20);
 line([0.158 0.218], [0.7 0.7], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');text(0.16,0.72,'* p = 0.01','FontName','sans serif','FontSize',20);
lh = legend(bl);
legnames = {'Post Sleep', 'Pre Sleep'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';
% hXLabel = xlabel('\fontsize{26} Time (s)');
% hYLabel = ylabel('\fontsize{26} auc');
set(gca,'FontName','sans serif','FontSize',14);
% set([ hXLabel, hYLabel, lh],'FontName','sans serif','FontSize',24);
ylim([0.5 0.75])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
%text(0.23,0.64,'p = 0.016','FontName','sans serif','FontSize',20);


%Wake - Wake
subplot(2,1,2);
bl=boundedline(ERP_preI3{i}.time,mean(decodeRI2_WG,2), std(decodeRI2_WG')./sqrt(17),ERP_preI3{i}.time,mean(decodeRI1_WG,2), std(decodeRI1_WG')./sqrt(17),'cmap', colors1, 'alpha','transparency', 0.2);
set(bl,'LineStyle','-','LineWidth',4) %set line width
line([0 0], [-250 250], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');
lh = legend(bl);
legnames = {'Post Wake', 'Pre Wake'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'northwest';

set(gca,'FontName','sans serif','FontSize',14);
% set([ hXLabel, hYLabel lh],'FontName','sans serif','FontSize',24);
ylim([0.5 0.75])
% Extra axis for boxing
haxes1 = gca; % handle to axes
haxes1_pos = get(haxes1,'Position'); % store position of first axes
set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

han=axes(f,'visible','off');
han.XLabel.Visible='on';
xlabel(han,sprintf('\n time (s)'), 'FontName','sans serif','FontSize',26);
han.YLabel.Visible='on';
ylabel(han,{'auc', ' '}, 'FontName','sans serif','FontSize',26);

%% I5

%  figure; colors1 = [0 0 1 ;1 0 0];
% % Sleep - Sleep
% bl=boundedline(ERP_preI5{i}.time,mean(decodeI45_SG,2), std(decodeI45_SG')./sqrt(16),ERP_preI5{i}.time,mean(decodeI45_WG,2), std(decodeI45_WG')./sqrt(17),'cmap', colors1, 'alpha','transparency', 0.2);
% set(bl,'LineStyle','-','LineWidth',4) %set line width
% line([0 0], [-250 250], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');
% % line([-0.054 0.06], [0.7 0.7], 'color','k','LineStyle','-', 'LineWidth',4, 'DisplayName', ''); text(-0.045,0.72,'*p = 0.01','FontName','sans serif','FontSize',20);
% % line([-0.18 -0.13], [0.7 0.7], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');text(-0.15,0.71,'+ p = 0.061','FontName','sans serif','FontSize',20);
% % line([-0.048 0.04], [0.7 0.7], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');text(0.01,0.71,'* p = 0.024','FontName','sans serif','FontSize',20);
% % line([0.25 0.284], [0.7 0.7], 'color','k','LineStyle','--', 'LineWidth',4, 'DisplayName', '');text(0.26,0.71,'+ p = 0.079','FontName','sans serif','FontSize',20);
% lh = legend(bl);
% legnames = {'Sleep', 'Wake'}; % here you can rename your conditions
% for i = 1:length(legnames)
%     str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
% end
% lh.String = str; lh.FontSize=16; lh.Box = 'off'; lh.Location = 'southwest';
% hXLabel = xlabel('\fontsize{26} Time (s)');
% hYLabel = ylabel('\fontsize{26} auc');
% set(gca,'FontName','sans serif','FontSize',14);
% set([ hXLabel, hYLabel, lh],'FontName','sans serif','FontSize',24);
% ylim([0.45 0.65])
% % Extra axis for boxing
% haxes1 = gca; % handle to axes
% haxes1_pos = get(haxes1,'Position'); % store position of first axes
% set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 
% %text(0.23,0.64,'p = 0.016','FontName','sans serif','FontSize',20);
% 
%% %% Within Subject Permutation - Searchlight
%%%%%%%%%%%%%%%%%
cfg=[];
cfg.statistic        ='ft_statfun_depsamplesT'; %between groups independant sampels statistics
cfg.method           ='montecarlo'; % Monte Carlo method for calculating the signif prob an estimate of the p-value under the permutation distribution.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
%  cfg.avgoverchan      = 'yes';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
cfg.uvar             = 1;
cfg.ivar             = 2;

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
cfg.design= design;  % design matrix

load('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\LAYOUT\ChanInfoFT.mat')
cfg.neighbours    = neighb;
% design the fieldtrip structures
cfg_g=[];
cfg_g.channel = 'C3';
cfg_g.keepindividual = 'yes';
pre_R2S=ft_timelockgrandaverage(cfg_g, ERP_preR2{SG});

Before = pre_R2S; Before.individual =  permute(decode_RI1S,[2,1]); Before.time = 1; Before.label=preR2.label; 
After  = pre_R2S; After.individual =  permute(decode_RI2S,[2,1]); After.time = 1; After.label=preR2.label;

Before_S=Before;Before_S.individual = Before.individual(SG,:,:); Before_S.dimord = 'subj_chan';
After_S=Before;After_S.individual = After.individual(SG,:,:);  After_S.dimord = 'subj_chan';

% run stats
[stat_SS]   = ft_timelockstatistics(cfg, After_S, Before_S);


%% NEW TOPO
load('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\LAYOUT\ChanInfoFT.mat')

temp          = resultRI_1S{1};
temp.perf     = mean(decodeRI1S_SG,2);
temp.perf_std = std(decodeRI1S_SG,0,2);

P   = temp.perf; 
idx = [1 2 3 4 11 5 6 7 8 9 10];
p1n = stat_SS.stat(idx);% new order of electrodes as data
h1n = stat_SS.negclusterslabelmat(idx);% new order of electrodes as data

% Create some dummy data
Data.label = lay.label; % Channel names
Data.dimord = 'chan_freq'; %dimensions of the data (channels x frequency)
% hier faket man quasi Fieldtrip vor, dass man mehrere Frequenzen hat,
% sonst macht er's nicht (d.h. wenn man einfach nur z.B. z-Werte plotten willl ?ber den Skalp)
Data.freq = 1:10; %arbitrary frequency vector
Data.zvalues = repelem(p1n,1,length(Data.freq)); 

% Do topoplot 
cfg = [];
cfg.layout = lay; %layout struct
cfg.parameter = 'zvalues';
cfg.colormap=colormap(brewermap(256, '*RdYlBu'));   
cfg.style = 'straight';
cfg.gridscale =800;
cfg.markersize =15;
cfg.highlightsize=60;
cfg.highlightcolor=[0 0 0];
cfg.highlightsymbol='.';
 cfg.zlim = [-3 3];
cfg.highlight = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';

cfg.highlightchannel = find(h1n == 1);
figure;ft_topoplotER(cfg,Data); title('Sleep (After - Before)'); 
set(gca,'FontSize',28,'TickLength',[.02 .02],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5); 

