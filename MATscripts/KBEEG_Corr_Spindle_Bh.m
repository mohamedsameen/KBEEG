% Spindle analysis useing the excel sheet produced from the siesta dominik algorithm
clear;

%fieldtrip
addpath(genpath('/home/b1044271/Toolboxes/obob_ownft'));
obob_init_ft;

cd '/home/b1044271/EEG_KB/Sleep/Spindles_detection/spindleFeatures/spindleResults'
path_R='/home/b1044271/EEG_KB/Sleep/Spindles_detection/spindleFeatures/spindleResults/';
dataRow_r   =struct2cell(dir(fullfile(path_R , '*.xls'))');
files_R =dataRow_r(1,:);


%% Groups 
AG = [1:10, 22, 23,24,25,31,32,33];
SG = [11:21, 26,27,28,29,30];

%% Performace Gain
PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);

PG_S = PG.Var1(SG); PG_W=PG.Var1(AG);

for i = 1:length(files_R)
    
    [~,sheet_name]=xlsfinfo([path_R files_R{i}]);


    for k=1:numel(sheet_name)
      data{k}=xlsread([path_R files_R{i}],sheet_name{k});
    end

    data2 = data(~cellfun('isempty',data));
    
    %% NREM
     NAmat = data2{1};
     NAmat_A = NAmat(1:2:end,:);
     NAmat_E = NAmat(2:2:end,:);

    %% N2
    N2mat = data2{2};
    N2mat_A = N2mat(1:2:end,:);
    N2mat_E = N2mat(2:2:end,:); 
    %% N3
    N3mat = data2{3};
     N3mat_A = N3mat(1:2:end,:);
     N3mat_E = N3mat(2:2:end,:);
     
    [r1(i),p1(i)] = corr(PG_S, NAmat_E(:,16)- NAmat_A(:,16), 'type','kendall');
    [r2(i),p2(i)] = corr(PG_S, N2mat_E(:,16)-N2mat_A(:,16), 'type','kendall');
    [r3(i),p3(i)] = corr(PG_S, N3mat_E(:,16)-N3mat_A(:,16), 'type','kendall');
    
%     Nall{} = NAmat_E(:,16)- NAmat_A(:,16)
%     writetable(table(NAmat_E(:,15)- NAmat_A(:,15), N2mat_E(:,15)- N2mat_A(:,15), N3mat_E(:,15)- N3mat_A(:,15), 'VariableNames',{'ALL','N2','N3'}), sprintf('%s_Fastspindles_Exp-Adpt.txt', files_R{i}(1:2)));
    
    
    
%     [r4(i),p4(i)] = corr(PG_S, NAmat_E(:,1), 'type','Kendall')
%     [r5,p5] = corr(PG_S, N2mat_E(:,15), 'type','Kendall')
%     [r6,p6] = corr(PG_S, N3mat_E(:,15), 'type','Kendall')
%     
% 
%     [r7,p7] =  corr(PG_S, NAmat_E(:,3), 'type','Kendall')
%     [r1,p1] = corr(PG_S, N2mat_E(:,3), 'type','Kendall')
%     [r1,p1] = corr(PG_S, N3mat_E(:,3), 'type','Kendall')
end



%% Topoplots

% Your data has to have the Fieldtrip format and you need to load channel
% layout info e.g.:
load('/home/b1044271/Toolboxes/LAYOUT/ChanInfoFT') % for my setup

idx = [4 6 5 1 3 2 9 11 10 7 8];

r1n = r1(idx); p1n = p1(idx);
r2n = r2(idx); p2n = p2(idx);
r3n = r3(idx); p3n = p3(idx);

% Create some dummy data
Data.label = lay.label; % Channel names
Data.dimord = 'chan_freq'; %dimensions of the data (channels x frequency)
% hier faket man quasi Fieldtrip vor, dass man mehrere Frequenzen hat,
% sonst macht er's nicht (d.h. wenn man einfach nur z.B. z-Werte plotten willl ?ber den Skalp)
Data.freq = 1:10; %arbitrary frequency vector
Data.zvalues = repelem(r1n',1,length(Data.freq)); 

Data2 = Data; Data2.zvalues = repelem(r2n',1,length(Data.freq)); 
Data3 = Data; Data3.zvalues = repelem(r3n',1,length(Data.freq)); 

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

[fdr,q] = mafdr(p1n)
% p1n2 = bonf_holm(p1n , 0.05);
cfg.highlightchannel = find(p1n*11 < 0.05);
figure; ft_topoplotER(cfg,Data); title('NREM'); 
set(gca,'FontSize',18,'TickLength',[.01 .01],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1.5); 

% p2n2 = bonf_holm(p2n , 0.05);
cfg.highlightchannel = find(p2n*11 < 0.05);
figure; ft_topoplotER(cfg,Data2); title('N2'); 
set(gca,'FontSize',18,'TickLength',[.01 .01],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1.5); 

% p3n2 = bonf_holm(p3n , 0.05);
cfg.highlightchannel = find(p3n*11 < 0.05);
figure; ft_topoplotER(cfg,Data3); title('N3'); 
set(gca,'FontSize',18,'TickLength',[.01 .01],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1.5); 


h = colorbar();
w = h.LineWidth;
h.LineWidth = 2;
ylabel(h, 't','FontSize', 22);
set(gca,'FontSize',18,'TickLength',[.01 .01],'XGrid','off','YGrid','off', 'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1.5); 



%% C3 C4 
%fieldtrip
% addpath(genpath('/home/b1044271/Toolboxes/obob_ownft'));
% obob_init_ft;

cd '/home/b1044271/EEG_KB/Sleep/Spindles_detection/spindleFeatures/spindleResults'
path_R='/home/b1044271/EEG_KB/Sleep/Spindles_detection/spindleFeatures/spindleResults/';
dataRow_r   =struct2cell(dir(fullfile(path_R , '*.xls'))');
files_R =dataRow_r(1,:);idx = [4 6 5 1 3 2 9 11 10 7 8];


for i = 1:2
    
    [~,sheet_name]=xlsfinfo([path_R files_R{i}]);


    for k=1:numel(sheet_name)
      data{k}=xlsread([path_R files_R{i}],sheet_name{k});
    end

    data2 = data(~cellfun('isempty',data));
    
    %% NREM
     NAmat= data2{1};
     NAmat_A(1:16,i) = NAmat(1:2:end,16);
     NAmat_E(1:16,i)  = NAmat(2:2:end,16);

    %% N2
    N2mat= data2{2};
    N2mat_A(1:16,i) = N2mat(1:2:end,16);
    N2mat_E(1:16,i) = N2mat(2:2:end,16); 
    %% N3
     N3mat= data2{3};
     N3mat_A(1:16,i) = N3mat(1:2:end,16);
     N3mat_E(1:16,i) = N3mat(2:2:end,16);
end   

SG = [11:21, 26,27,28,29,30];

%Performace Gain
PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);

PG_S = PG.Var1(SG); 
    [r1,p1] = corr(PG_S, mean(NAmat_E,2)- mean(NAmat_A,2), 'type','kendall')
    [r2,p2] = corr(PG_S, mean(N2mat_E,2)-mean(N2mat_A,2), 'type','kendall')
    [r3,p3] = corr(PG_S, mean(N3mat_E,2)-mean(N3mat_A,2), 'type','kendall')
    
    
 writetable(table(mean(NAmat_E,2)- mean(NAmat_A,2), mean(N2mat_E,2)-mean(N2mat_A,2), mean(N3mat_E,2)-mean(N3mat_A,2), ...
     'VariableNames',{'ALL','N2','N3'}), 'Fastspindles_Exp-Adpt_density_central.txt');


 writetable(table(mean(NAmat_E,2), mean(N2mat_E,2), mean(N3mat_E,2), ...
     'VariableNames',{'ALL','N2','N3'}), 'Fastspindles_ExpOnly_density_central.txt');
%% CLUSTER_BASED Permutation test for CORRELATIONS

clear;

%fieldtrip
addpath(genpath('/home/b1044271/Toolboxes/obob_ownft'));
obob_init_ft;

% path for data
cd '/home/b1044271/EEG_KB/Sleep/Spindles_detection/spindleFeatures/spindleResults'
path_R='/home/b1044271/EEG_KB/Sleep/Spindles_detection/spindleFeatures/spindleResults/';
dataRow_r   =struct2cell(dir(fullfile(path_R , '*.xls'))');
files_R =dataRow_r(1,:);

% Path for tools
addpath(genpath('/home/b1044271/Toolboxes/swtest')); % normality test
addpath(genpath('/home/b1044271/Toolboxes/cluster_corr')); % correlation permutation function

% neighbours electrodes
load('/home/b1044271/Toolboxes/LAYOUT/ChanInfoFT') % for my setup

% Performace Gain (Behavior)
PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);
SG = [11:21, 26,27,28,29,30];PG_S = PG.Var1(SG); 

% load files
for i = 1:length(files_R)
    
    [~,sheet_name]=xlsfinfo([path_R files_R{i}]);


    for k=1:numel(sheet_name)
      data{k}=xlsread([path_R files_R{i}],sheet_name{k}); % read all the sheets in the excel file
    end

    data2 = data(~cellfun('isempty',data));
    
     %% NREM
     NAmat= data2{1};
     NAmat_A(1:16,i) = NAmat(1:2:end,16); % Column 16 is the density of fast spindles
     NAmat_E(1:16,i)  = NAmat(2:2:end,16);

    %% N2
    N2mat= data2{2};
    N2mat_A(1:16,i) = N2mat(1:2:end,16);
    N2mat_E(1:16,i) = N2mat(2:2:end,16); 
    %% N3
     N3mat= data2{3};
     N3mat_A(1:16,i) = N3mat(1:2:end,16);
     N3mat_E(1:16,i) = N3mat(2:2:end,16);
end



% load the channel data; needs fieldtrip compatible channel layout
% (ChanInfoFT.mat)
load('/home/b1044271/Toolboxes/LAYOUT/ChanInfoFT') % for my setup

PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);
SG = [11:21, 26,27,28,29,30];PG_S = PG.Var1(SG); 

idx = [4 6 5 1 3 2 9 11 10 7 8]; % the correct order of the spindle results files 

% N2
corr_data = N2mat_E(:,idx)- N2mat_A(:,idx); % put the spindle results in the correct order as the layout of the 11 eelectrodes
stat_fs = clustertest_corr(corr_data, PG_S, 0.05, 0.025, 1);
plot_topo_mo(stat_fs.rho, stat_fs.posclusterslabelmat, 'N2')

%N3
corr_data = N3mat_E(:,idx)- N3mat_A(:,idx); % put the spindle results in the correct order as the layout of the 11 eelectrodes
stat_fs2 = clustertest_corr(corr_data, PG_S, 0.05, 0.025, 1);
plot_topo_mo(stat_fs.rho, stat_fs.posclusterslabelmat, 'N3')

%NREM
corr_data = NAmat_E(:,idx)- NAmat_A(:,idx); % put the spindle results in the correct order as the layout of the 11 eelectrodes
stat_fs3 = clustertest_corr(corr_data, PG_S, 0.05, 0.025, 1);
plot_topo_mo(stat_fs.rho, stat_fs.posclusterslabelmat, 'NREM')


%% COHEN D

%GRAND AVERAGE (optional)

grandavgFIC_sel.powspctrm = corr_data;


grandavgFC_sel.powspctrm = repmat(PG_S,1,11);

% COHEN'S D
x1 = nan(16,1);
x2 = nan(16,1);

for i=1:16

  %construct a 3-dimensional Boolean array to select the data from this participant
  sel3d = false(size(grandavgFIC_sel.powspctrm));
  sel3d(i,:,:) = stat_fs .posclusterslabelmat==1;

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