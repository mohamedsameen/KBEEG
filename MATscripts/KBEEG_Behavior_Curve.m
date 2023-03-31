% Behavioral analysis KB EEG experiment
% metric --> performance = accuracy (correct letters per word) / time (s)
% for each subject one value for each sessions/folder (INV1-5 - R1-2)
%look at the results by letter file in each folder

clear all;
%fieldtrip
addpath(genpath('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\obob_ownft'));
obob_init_ft;

AG = [1:10, 22, 23,24,25,31,32,33];
SG = [11:21, 26,27,28,29,30];
addpath(genpath('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Toolboxes\kakearney-boundedline-pkg-50f7e4b'))
%% load files
cd C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\BehavioralData;
path = 'C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\BehavioralData\';

% get all the VP folders in the directory
vps=dir(path);
idx=find(~cellfun(@isempty,regexp({vps(:).name}, 'VP'))); 
VPs={vps(idx).name};

opts = detectImportOptions('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\BehavioralData\VP01\VP01INV1_RESULTS_BY_WORD.txt');


pr_INV1=nan(1000,1000);
pr_INV2=nan(1000,1000);
pr_INV3=nan(1000,1000);
pr_INV4=nan(1000,1000);
pr_INV5=nan(1000,1000);
pr_REG1=nan(1000,1000);
pr_REG2=nan(1000,1000);

for subj = 1:length(VPs)
    
    INV1 = readtable([path VPs{subj} '\' VPs{subj} 'INV1_RESULTS_BY_WORD.txt'], opts);
    INV2 = readtable([path VPs{subj} '\' VPs{subj} 'INV2_RESULTS_BY_WORD.txt'], opts);
    INV3 = readtable([path VPs{subj} '\' VPs{subj} 'INV3_RESULTS_BY_WORD.txt'], opts);
    INV4 = readtable([path VPs{subj} '\' VPs{subj} 'INV4_RESULTS_BY_WORD.txt'], opts);
    INV5 = readtable([path VPs{subj} '\' VPs{subj} 'INV5_RESULTS_BY_WORD.txt'], opts);
    
    REG1 = readtable([path VPs{subj} '\' VPs{subj} 'REG1_RESULTS_BY_WORD.txt'], opts);
    REG2 = readtable([path VPs{subj} '\' VPs{subj} 'REG2_RESULTS_BY_WORD.txt'], opts);
    
    % Length of session
    INV1_l=length(table2array(INV1(:,8)));
    INV2_l=length(table2array(INV2(:,8)));
    INV3_l=length(table2array(INV3(:,8)));
    INV4_l=length(table2array(INV4(:,8)));
    INV5_l=length(table2array(INV5(:,8)));
    
    REG1_l=length(table2array(REG1(:,8)));
    REG2_l=length(table2array(REG2(:,8)));
    
    %% Accuracy
    acc_INV1(1:INV1_l)  =  table2array(INV1(:,8))./5;
    acc_INV2(1:INV2_l)  =  table2array(INV2(:,8))./5;
    acc_INV3(1:INV3_l)  =  table2array(INV3(:,8))./5;
    acc_INV4(1:INV4_l)  =  table2array(INV4(:,8))./5;
    acc_INV5(1:INV5_l)  =  table2array(INV5(:,8))./5;
    
    acc_REG1(1:REG1_l)  =  table2array(REG1(:,8))./5;
    acc_REG2(1:REG2_l)  =  table2array(REG2(:,8))./5;
    
    %% speed
    sp_INV1(1:INV1_l)=table2array(INV1(:,10))./1000;
    sp_INV2(1:INV2_l)=table2array(INV2(:,10))./1000;
    sp_INV3(1:INV3_l)=table2array(INV3(:,10))./1000;
    sp_INV4(1:INV4_l)=table2array(INV4(:,10))./1000;
    sp_INV5(1:INV5_l)=table2array(INV5(:,10))./1000;
    
    sp_REG1(1:REG1_l)=table2array(REG1(:,10))./1000;
    sp_REG2(1:REG2_l)=table2array(REG2(:,10))./1000;
    
    
    %% Performance (accuracy/time)
    pr_INV1(1:INV1_l,subj)=acc_INV1(:)./sp_INV1(:); 
    pr_INV2(1:INV2_l,subj)=acc_INV2(:)./sp_INV2(:);
    pr_INV3(1:INV3_l,subj)=acc_INV3(:)./sp_INV3(:);
    pr_INV4(1:INV4_l,subj)=acc_INV4(:)./sp_INV4(:);
    pr_INV5(1:INV5_l,subj)=acc_INV5(:)./sp_INV5(:);
    
    pr_REG1(1:REG1_l,subj)=acc_REG1(:)./sp_REG1(:);
    pr_REG2(1:REG2_l,subj)=acc_REG2(:)./sp_REG2(:);
   
    clearvars acc_INV1 acc_INV2 acc_INV3 acc_INV4 acc_INV5 sp_INV1 sp_INV2 sp_INV3 sp_INV4 sp_INV5 sp_REG1 sp_REG2 acc_REG1 acc_REG2
end

for i = 1:33
    
    B1 = ~isnan(pr_INV1(:,i));
    B2 = ~isnan(pr_INV2(:,i));
    B3 = ~isnan(pr_INV3(:,i));
    B4 = ~isnan(pr_INV4(:,i));
    B5 = ~isnan(pr_INV5(:,i));
    B6 = ~isnan(pr_REG1(:,i));
    B7 = ~isnan(pr_REG2(:,i));

    Indices1(i) = arrayfun(@(x) find(B1(:, x), 1, 'last'), 1:size(pr_INV1(:,i), 2));
    Indices2(i) = arrayfun(@(x) find(B2(:, x), 1, 'last'), 1:size(pr_INV2(:,i), 2));
    Indices3(i) = arrayfun(@(x) find(B3(:, x), 1, 'last'), 1:size(pr_INV3(:,i), 2));
    Indices4(i) = arrayfun(@(x) find(B4(:, x), 1, 'last'), 1:size(pr_INV4(:,i), 2));
    Indices5(i) = arrayfun(@(x) find(B5(:, x), 1, 'last'), 1:size(pr_INV5(:,i), 2));
    Indices6(i) = arrayfun(@(x) find(B6(:, x), 1, 'last'), 1:size(pr_REG1(:,i), 2));
    Indices7(i) = arrayfun(@(x) find(B7(:, x), 1, 'last'), 1:size(pr_REG2(:,i), 2));

end

X1s = pr_INV1(1:min(Indices1),SG); X1w = pr_INV1(1:min(Indices1),AG);
X2s = pr_INV2(1:min(Indices2),SG); X2w = pr_INV2(1:min(Indices2),AG);
X3s = pr_INV3(1:min(Indices3),SG); X3w = pr_INV3(1:min(Indices3),AG);
X4s = pr_INV4(1:min(Indices4),SG); X4w = pr_INV4(1:min(Indices4),AG);
X5s = pr_INV5(1:min(Indices5),SG); X5w = pr_INV5(1:min(Indices5),AG);
X6s = pr_REG1(1:min(Indices6),SG); X6w = pr_REG1(1:min(Indices6),AG);
X7s = pr_REG2(1:min(Indices7),SG); X7w = pr_REG2(1:min(Indices7),AG);

X34s = [X3s;X4s]; X34w = [X3w;X4w];

colors1 = [1 0 0; 0 0 1];
IF = figure;

subplot(2,3,1)
yI1 = smooth(1:size(X1w,1),nanmean(X1w,2), 0.1, 'loess'); % smoothing
yI12 = smooth(1:size(X1w,1),nanmean(X1s,2), 0.1, 'loess'); % smoothing
yI1s = smooth(1:size(X1w,1),nanstd(X1w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI12s = smooth(1:size(X1w,1),nanstd(X1s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing
b1=boundedline(1:size(X1w,1),yI1, yI1s,1:size(X1w,1),yI12,yI12s,...
    'cmap', colors1, 'alpha','transparency', 0.2);
xlim([0 size(X1w,1)]); ylim([0 0.4]);set(b1,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

subplot(2,3,2)
yI2 = smooth(1:size(X2w,1),nanmean(X2w,2), 0.1, 'loess'); % smoothing
yI22 = smooth(1:size(X2w,1),nanmean(X2s,2), 0.1, 'loess'); % smoothing
yI2s = smooth(1:size(X2w,1),nanstd(X2w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI22s = smooth(1:size(X2w,1),nanstd(X2s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing
b2=boundedline(1:size(X2w,1),yI2, yI2s,1:size(X2w,1),yI22,yI22s,...
    'cmap', colors1, 'alpha','transparency', 0.2);
xlim([0 size(X2w,1)]); ylim([0 0.4]);set(b2,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

subplot(2,3,3)
yI3 = smooth(1:size(X3w,1),nanmean(X3w,2), 0.1, 'loess'); % smoothing
yI32 = smooth(1:size(X3w,1),nanmean(X3s,2), 0.1, 'loess'); % smoothing
yI3s = smooth(1:size(X3w,1),nanstd(X3w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI32s = smooth(1:size(X3w,1),nanstd(X3s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing
b3=boundedline(1:size(X3w,1),yI3, yI3s,1:size(X3w,1),yI32,yI32s,...
    'cmap', colors1, 'alpha','transparency', 0.2);
xlim([0 size(X3w,1)]); ylim([0 0.4]);set(b3,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

subplot(2,3,4)
yI4 = smooth(1:size(X4w,1),nanmean(X4w,2), 0.1, 'loess'); % smoothing
yI42 = smooth(1:size(X4w,1),nanmean(X4s,2), 0.1, 'loess'); % smoothing
yI4s = smooth(1:size(X4w,1),nanstd(X4w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI42s = smooth(1:size(X4w,1),nanstd(X4s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing
b4=boundedline(1:size(X4w,1),yI4, yI4s,1:size(X4w,1),yI42,yI42s,...
    'cmap', colors1, 'alpha','transparency', 0.2);
xlim([0 size(X4w,1)]); ylim([0 0.4]);set(b4,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

subplot(2,3,5)
yI5 = smooth(1:size(X5w,1),nanmean(X5w,2), 0.1, 'loess'); % smoothing
yI52 = smooth(1:size(X5w,1),nanmean(X5s,2), 0.1, 'loess'); % smoothing
yI5s = smooth(1:size(X5w,1),nanstd(X5w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI52s = smooth(1:size(X5w,1),nanstd(X5s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing
b5=boundedline(1:size(X5w,1),yI5, yI5s,1:size(X5w,1),yI52,yI52s,...
    'cmap', colors1, 'alpha','transparency', 0.2);
xlim([0 size(X5w,1)]); ylim([0 0.4]);set(b5,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

AxesH    = findobj(IF, 'Type', 'Axes');
YLabelHC = get(AxesH, 'YLabel');
YLabelH  = [YLabelHC{:}];
set(YLabelH, 'String', 'performance')
XLabelHC = get(AxesH, 'XLabel');
XLabelH  = [XLabelHC{:}];
set(XLabelH, 'String', 'trials (n)')
lh = legend([b1; b2; b3 ;b4 ;b5]);legnames = {'Wake', 'Sleep'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off'; 
%%%%%%%%%%%%%%%%%%%

%% REG
colors1 = [1 0 0; 0 0 1];
IF = figure;

subplot(2,1,1)
yI1 = smooth(1:size(X6w,1),nanmean(X6w,2), 0.1, 'loess'); % smoothing
yI12 = smooth(1:size(X6w,1),nanmean(X6s,2), 0.1, 'loess'); % smoothing
yI1s = smooth(1:size(X6w,1),nanstd(X6w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI12s = smooth(1:size(X6w,1),nanstd(X6s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing
b1=boundedline(1:size(X6w,1),yI1, yI1s,1:size(X6w,1),yI12,yI12s,...
    'cmap', colors1, 'alpha','transparency', 0.2);
xlim([1 size(X6w,1)]); ylim([0.4 1.1]);set(b1,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

subplot(2,1,2)
yI2 = smooth(1:size(X7w,1),nanmean(X7w,2), 0.1, 'loess'); % smoothing
yI22 = smooth(1:size(X7w,1),nanmean(X7s,2), 0.1, 'loess'); % smoothing
yI2s = smooth(1:size(X7w,1),nanstd(X7w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI22s = smooth(1:size(X7w,1),nanstd(X7s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing
b2=boundedline(1:size(X7w,1),yI2, yI2s,1:size(X7w,1),yI22,yI22s,...
    'cmap', colors1, 'alpha','transparency', 0.2);
xlim([1 size(X7w,1)]); ylim([0.4 1.1]);set(b2,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

AxesH    = findobj(IF, 'Type', 'Axes');
YLabelHC = get(AxesH, 'YLabel');
YLabelH  = [YLabelHC{:}];
set(YLabelH, 'String', 'performance')
XLabelHC = get(AxesH, 'XLabel');
XLabelH  = [XLabelHC{:}];
set(XLabelH, 'String', 'trials (n)')
lh = legend([b1; b2]);legnames = {'Wake', 'Sleep'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off';
%% Permutation

% design the fieldtrip structures
load('C:\Users\b1044271\Desktop\erp_FV.mat')
cfg=[];
cfg_g.channel = 'E21';
cfg_g.keepindividual = 'yes';
pre_R2S=ft_timelockgrandaverage(cfg_g, erp_FV{:});
pre_R2S.time = 1:1:34;

Sleep_S = pre_R2S; Sleep_S.individual =  permute(X4s,[2,3,1]); Sleep_S.time = 1:length(X4s);
Wake_S = pre_R2S; Wake_S.individual =  permute(X4w,[2,3,1]); Wake_S.time = 1:length(X4s);

Sleep_S2 = pre_R2S; Sleep_S2.individual =  permute(X5s,[2,3,1]); Sleep_S2.time = 1:length(X5s);
Wake_S2 = pre_R2S; Wake_S2.individual =  permute(X5w,[2,3,1]); Wake_S2.time = 1:length(X5s);
%permutation
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

[stat_S]   = ft_timelockstatistics(cfg, Sleep_S, Wake_S);
[stat_S2]   = ft_timelockstatistics(cfg, Sleep_S2, Wake_S2);

%% COHEN D
%%%%%%%%%%
cfg=[];
cfg.statistic        ='ft_statfun_indepsamplesT'; %between groups independant sampels statistics
cfg.method           ='montecarlo'; % Monte Carlo method for calculating the signif prob an estimate of the p-value under the permutation distribution.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.avgoverchan      = 'yes';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 100000;
cfg.correcttail = 'prob';
cfg.ivar             = 1;

% Design Matrix for T-Test (2 Conditions)
design = [ones(1,16), ones(1,17)*2];
cfg.design= design;  

inference = ft_timelockstatistics(cfg, Sleep_S, Wake_S);

% GRAND AVERAGE:
cfg = [];
cfg.latency     = [-0.5 2];
grandavgFIC_sel = Sleep_S;
grandavgFC_sel  = Wake_S;

x1 = nan(17,1);
x2 = nan(17,1);

for i=1:17

  % construct a 3-dimensional Boolean array to select the data from this participant
  sel3d = false(size(grandavgFIC_sel.individual));
  sel3d(i,:,:) = inference.posclusterslabelmat==1;

  % select the FIC data in the cluster for this participant, represent it as a vector
  tmp = grandavgFIC_sel.individual(sel3d(:));
  % compute the average over the cluster
  x1(i) = mean(tmp);

  % select the FC data in the cluster for this participant, represent it as a vector
  tmp = grandavgFC_sel.individual(sel3d(:));
  % compute the average over the cluster
  x2(i) = mean(tmp);
end

n1 = length(x1);
n2 = length(x2);


cohensd = mean(x1-x2) ./ std(x1-x2);

% writetable(table(X4s), 'Sleep_BH_trial.txt', 'WriteVariableNames',0)
% writetable(table(X4w), 'Wake_BH_trial.txt', 'WriteVariableNames',0)

% save('BH_trial_S1.mat','X1s')
% save('BH_trial_W1.mat','X1w')
% save('BH_trial_S2.mat','X2s')
% save('BH_trial_W2.mat','X2w')
% save('BH_trial_S3.mat','X3s')
% save('BH_trial_W3.mat','X3w')
% save('BH_trial_S5.mat','X5s')
% save('BH_trial_W5.mat','X5w')
% save('BH_trial_S.mat','X4s')
% save('BH_trial_W.mat','X4w')


save('BH_trial_S34.mat','X34s')
save('BH_trial_W34.mat','X34w')

%% LISTS IN REG

X = readtable('C:\Users\b1044271\Desktop\Sodium\Projects\Keyboard\EEG\Results\2022\my_data.csv');
X.Var1=[];
AG2 = [1:7,9:10, 22, 23,24,25,31,32,33];

XS = table2array(X(:,SG));
XW = table2array(X(:,AG));
XW1= str2double(table2array(X(1:12,8)));
XW(1:12,17)=XW1(1:12);

colors1 = [1 0 0; 0 0 1];
IF = figure;

yI2 = smooth(1:size(XW,1),nanmean(XW,2), 0.1, 'loess'); % smoothing
yI22 = smooth(1:size(XS,1),nanmean(XS,2), 0.1, 'loess'); % smoothing
yI2s = smooth(1:size(XW,1),nanstd(XW,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI22s = smooth(1:size(XS,1),nanstd(XS,[],2)./sqrt(16), 0.1, 'loess'); % smoothing
b2=boundedline(1:size(XW,1),yI2, yI2s,1:size(XS,1),yI22,yI22s,...
    'cmap', colors1, 'alpha','transparency', 0.2);
xlim([1 size(XW,1)]); ylim([0.4 1.1]);set(b2,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

AxesH    = findobj(IF, 'Type', 'Axes');
YLabelHC = get(AxesH, 'YLabel');
YLabelH  = [YLabelHC{:}];
set(YLabelH, 'String', 'performance')
XLabelHC = get(AxesH, 'XLabel');
XLabelH  = [XLabelHC{:}];
set(XLabelH, 'String', 'trials (n)')
lh = legend([b1; b2]);legnames = {'Wake', 'Sleep'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off';

%% Permutation List 1

% design the fieldtrip structures
load('C:\Users\b1044271\Desktop\erp_FV.mat')
cfg=[];
cfg_g.channel = 'E21';
cfg_g.keepindividual = 'yes';
pre_R2S=ft_timelockgrandaverage(cfg_g, erp_FV{:});
pre_R2S.time = 1:1:12;

Sleep_Sr = pre_R2S; Sleep_Sr.individual =  permute(XS,[2,3,1]); Sleep_Sr.time = 1:size(XS,1);
Wake_Sr = pre_R2S; Wake_Sr.individual =  permute(XW,[2,3,1]); Wake_Sr.time = 1:size(XW,1);


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

[stat_S]   = ft_timelockstatistics(cfg, Sleep_Sr, Wake_Sr);


disp(cohensd)

% save('R2L1_S.mat','XS');
% save('R2L1_W.mat','XW');

%% PLOT I3 and I4 together
colors1 = [1 0 0; 0 0 1];
IF = figure;

yI3 = smooth(1:size(X3w,1),nanmean(X3w,2), 0.1, 'loess'); % smoothing
yI32 = smooth(1:size(X3w,1),nanmean(X3s,2), 0.1, 'loess'); % smoothing
yI3s = smooth(1:size(X3w,1),nanstd(X3w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI32s = smooth(1:size(X3w,1),nanstd(X3s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing

yI4 = smooth(1:size(X4w,1),nanmean(X4w,2), 0.1, 'loess'); % smoothing
yI42 = smooth(1:size(X4w,1),nanmean(X4s,2), 0.1, 'loess'); % smoothing
yI4s = smooth(1:size(X4w,1),nanstd(X4w,[],2)./sqrt(17), 0.1, 'loess'); % smoothing
yI42s = smooth(1:size(X4w,1),nanstd(X4s,[],2)./sqrt(16), 0.1, 'loess'); % smoothing

yI34 = [yI3 ;yI4];
yI342 = [yI32;yI42];
yI34s = [yI3s ;yI4s];
yI342s = [yI32s;yI42s];

bX=boundedline(1:size(yI34,1),yI34, yI34s,1:size(yI34,1),yI342,yI342s,...
    'cmap', colors1, 'alpha','transparency', 0.2);

xlim([0 size(yI34,1)]); ylim([0 0.4]);set(bX,'LineStyle','-','LineWidth',4) %set line width
set(gca,'FontName','sans serif','FontSize',16); set(gca,'Box','off','TickLength',[.0 .0],...
    'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',2.5, 'FontSize', 20); 

line([73.5 73.5], [0. 0.4], 'color','k','LineStyle','--', 'LineWidth',3, 'DisplayName', '');

AxesH    = findobj(IF, 'Type', 'Axes');
YLabelHC = get(AxesH, 'YLabel');
YLabelH  = [YLabelHC{:}];
set(YLabelH, 'String', 'performance')
XLabelHC = get(AxesH, 'XLabel');
XLabelH  = [XLabelHC{:}];
set(XLabelH, 'String', 'trials (n)')
lh = legend(bX);legnames = {'Wake', 'Sleep'}; % here you can rename your conditions
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors1(i, 1), colors1(i, 2), colors1(i, 3), legnames{i})];
end
lh.String = str; lh.FontSize=16; lh.Box = 'off';

