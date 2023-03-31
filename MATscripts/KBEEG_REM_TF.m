%fieldtrip
clear;
addpath(genpath('/home/b1044271/Toolboxes/obob_ownft'));
obob_init_ft;

%EEGLab
addpath(genpath('/home/b1044271/Toolboxes/eeglab14_1_1b'));
eeglab; close gcf;


%% load EEG

path_E='/home/b1044271/EEG_KB/Sleep/Sleep_LPSG/';
dataRow_r   =struct2cell(dir(fullfile(path_E , '*.set'))');
files_E =dataRow_r(1,:);

path_A='/home/b1044271/EEG_KB/Sleep/Sleep_ADAPT/';
dataRow_r   =struct2cell(dir(fullfile(path_A , '*.set'))');
files_A =dataRow_r(1,:);

%% LOAD STAGES

path_SE = '/home/b1044271/EEG_KB/Sleep/Stages/Markers/Exp/';
dataRow_r   =struct2cell(dir(fullfile(path_SE , '*.txt'))');
St_E =dataRow_r(1,:);

path_SA = '/home/b1044271/EEG_KB/Sleep/Stages/Markers/Adapt/';
dataRow_r   =struct2cell(dir(fullfile(path_SA , '*.txt'))');
St_A =dataRow_r(1,:);

%% LOAD REM INFO

% ADAPTATION P/T
path_AREM = '/home/b1044271/EEG_KB/Sleep/REM/Markers/A_REM/';
dataRow_r   =struct2cell(dir(fullfile(path_AREM, '*.txt'))');
rem_A =dataRow_r(1,:);

% EXP PHASIC/TONIC
path_EREM = '/home/b1044271/EEG_KB/Sleep/REM/Markers/E_REM/';
dataRow_r   =struct2cell(dir(fullfile(path_EREM, '*.txt'))');
rem_E =dataRow_r(1,:);

% ADAPTATION PHASIC/TONIC
path_AEM = '/home/b1044271/EEG_KB/Sleep/REM/Markers/A_EM/';
dataRow_r   =struct2cell(dir(fullfile(path_AEM, '*.txt'))');
em_A =dataRow_r(1,:);

% EXP EMs
path_EEM = '/home/b1044271/EEG_KB/Sleep/REM/Markers/E_EM/';
dataRow_r   =struct2cell(dir(fullfile(path_EEM, '*.txt'))');
em_E =dataRow_r(1,:);

%% Groups 
AG = [1:10, 22, 23,24,25,31,32,33];
SG = [11:21, 26,27,28,29,30];

%% START
for i = 1:length(files_E)

    EEGE    = pop_loadset([path_E files_E{i}]);
    STAGE_E = [path_SE St_E{i}];
%     REM_E   = [path_EREM rem_E{i}];
    EM_E    = [path_EEM em_E{i}];
    
    EEGA=pop_loadset([path_A files_A{i}]);
    STAGE_A = [path_SA St_A{i}];
% REM_A   = [path_AREM rem_A{i}];
    EM_A    = [path_AEM em_A{i}];
    
    %% Load markers
    eventstruct=importevent(STAGE_E,EEGE.event,EEGE.srate,'fields',{'type','latency'},'timeunit',NaN);EEGE.event=eventstruct;
    eventstruct=importevent(STAGE_A,EEGA.event,EEGA.srate,'fields',{'type','latency'},'timeunit',NaN);EEGA.event=eventstruct;
    
%     eventstruct=importevent(REM_E,EEGE.event,EEGE.srate,'fields',{'type','latency'},'timeunit',NaN);EEGE.event=eventstruct;
%     eventstruct=importevent(REM_A,EEGA.event,EEGA.srate,'fields',{'type','latency'},'timeunit',NaN);EEGA.event=eventstruct;   
    
    eventstruct=importevent(EM_E,EEGE.event,EEGE.srate,'fields',{'type','latency'},'timeunit',NaN);EEGE.event=eventstruct;
    eventstruct=importevent(EM_A,EEGA.event,EEGA.srate,'fields',{'type','latency'},'timeunit',NaN);EEGA.event=eventstruct;

    %% Prepro
    % delete mastoids, eye and muscles
    EEGE = pop_select(EEGE, 'channel', 1:10); 
    EEGA = pop_select(EEGA, 'channel', 1:10); 
    
%     EEGE = pop_resample(EEGE,125);
%     EEGA = pop_resample(EEGA,125);

    %add the REF channel back and perform the average reference  
    EEGE                                   = pop_eegfiltnew(EEGE, 0.1, 60); %pop_eegfiltnew(EEG1, 13, 30);
%     EEGE = pop_cleanline(EEGE, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEGE.nbchan                            = EEGE.nbchan+1;
    EEGE.data(end+1,:)                     = zeros(1, EEGE.pnts);
    EEGE.chanlocs(1,EEGE.nbchan).labels    = 'CZ';
    EEGE                                   = pop_chanedit(EEGE, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEGE                                   = eeg_checkset( EEGE);
    EEGE                                   = pop_reref( EEGE, []);  % normal eeglab reref    
    
    %add the REF channel back and perform the average reference 
    EEGA                                   = pop_eegfiltnew(EEGA, 0.1, 60); %pop_eegfiltnew(EEG1, 13, 30);
%     EEGA = pop_cleanline(EEGA, 'bandwidth',2,'chanlist',[1:10] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEGA.nbchan                            = EEGA.nbchan+1;
    EEGA.data(end+1,:)                     = zeros(1, EEGA.pnts);
    EEGA.chanlocs(1,EEGA.nbchan).labels    = 'CZ';
    EEGA                                   = pop_chanedit(EEGA, 'lookup','/home/b1044271/Toolboxes/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','setref',{'1:11' 'CZ'});
    EEGA                                   = eeg_checkset( EEGA);
    EEGA                                   = pop_reref( EEGA, []);  % normal eeglab reref    
    
    
     %% segment night in half
%     EMAX = EEGE.xmax./2;
%     EEGE2=pop_select(EEGE,'time',[EMAX EEGE.xmax]);
% 
%     AMAX = EEGA.xmax./2;
%     EEGA2=pop_select(EEGA,'time',[AMAX EEGA.xmax]);

    %% Epoching
    
    % PHASIC & TONIC
%     pREM_E=pop_epoch(EEGE,{'51'},[-30 0]); pRs(i,1) = length(pREM_E.epoch);
%     pREM_A=pop_epoch(EEGA,{'51'},[-30 0]); pRs(i,2) = length(pREM_A.epoch);
%     
%     tREM_E=pop_epoch(EEGE,{'52'},[-30 0]); tRs(i,1) = length(tREM_E.epoch);
%     tREM_A=pop_epoch(EEGA,{'52'},[-30 0]); tRs(i,2) = length(tREM_A.epoch);
    
    % REMs & SEMs
    R_E=pop_epoch(EEGE,{'stage_5'},[-30 0]); Rstages(i,1) = length(R_E.epoch);
    R_A=pop_epoch(EEGA,{'stage_5'},[-30 0]); Rstages(i,2) = length(R_A.epoch);
    
%     REMs_E=pop_epoch(R_E,{'R'},[-2 2]); RRs(i,1) = length(REMs_E.epoch);
%     REMs_A=pop_epoch(R_A,{'R'},[-2 2]); RRs(i,2) = length(REMs_A.epoch);
%     
%     NEMs_E=pop_epoch(R_E,{'N'},[-2 2]); NRs(i,1) = length(NEMs_E.epoch);
%     NEMs_A=pop_epoch(R_A,{'N'},[-2 2]); NRs(i,2) = length(NEMs_A.epoch);
%     
    %% Equalize trials
%     minP = min(length(pREM_E.epoch), length(pREM_A.epoch));
%         mint = min(length(tREM_E.epoch), length(tREM_A.epoch));
%             minR = min(length(REMs_E.epoch), length(REMs_A.epoch));
%             minN = min(length(NEMs_E.epoch), length(NEMs_A.epoch));
     minS = min(length( R_E.epoch), length(R_A.epoch));

%     pREM_E2=pop_select(pREM_E, 'trial', datasample(1:length(pREM_E.epoch),minP,'Replace',false));
%     pREM_A2=pop_select(pREM_A, 'trial', datasample(1:length(pREM_A.epoch),minP,'Replace',false));
%     tREM_E2=pop_select(tREM_E, 'trial', datasample(1:length(tREM_E.epoch),mint,'Replace',false));
%     tREM_A2=pop_select(tREM_A, 'trial', datasample(1:length(tREM_A.epoch),mint,'Replace',false));
%     REMs_E2=pop_select(REMs_E, 'trial', datasample(1:length(REMs_E.epoch),minR,'Replace',false));
%     REMs_A2=pop_select(REMs_A, 'trial', datasample(1:length(REMs_A.epoch),minR,'Replace',false));
%     NEMs_E2=pop_select(NEMs_E, 'trial', datasample(1:length(NEMs_E.epoch),minN,'Replace',false));
%     NEMs_A2=pop_select(NEMs_A, 'trial', datasample(1:length(NEMs_A.epoch),minN,'Replace',false));
    Rs_E2=pop_select(R_E, 'trial', datasample(1:length(R_E.epoch),minS,'Replace',false));
    Rs_A2=pop_select(R_A, 'trial', datasample(1:length(R_A.epoch),minS,'Replace',false));
  
    %% TF analysis
    
    EEG1f=eeglab2fieldtrip(Rs_E2,'preprocessing','none');
    EEG2f=eeglab2fieldtrip(Rs_A2,'preprocessing','none');    
%     
%     EEG3f=eeglab2fieldtrip(tREM_E2,'preprocessing','none');
%     EEG4f=eeglab2fieldtrip(tREM_A2,'preprocessing','none');    
    
%     EEG5f=eeglab2fieldtrip(REMs_E2,'preprocessing','none');
%     EEG6f=eeglab2fieldtrip(REMs_A2,'preprocessing','none');
% %     
%     EEG7f=eeglab2fieldtrip(NEMs_E2,'preprocessing','none');
%     EEG8f=eeglab2fieldtrip(NEMs_A2,'preprocessing','none');

    cfg                    = [];
    cfg.foi               = 1:0.5:90;
    cfg.pad             = 'nextpow2';
    cfg.method       = 'irasa';
    
    cfg.output        = 'fractal';
    EEG1a{i}=ft_freqanalysis(cfg, EEG1f);
    EEG2a{i}=ft_freqanalysis(cfg, EEG2f);
%    EEG5a{i}=ft_freqanalysis(cfg, EEG5f);
%    EEG6a{i}=ft_freqanalysis(cfg, EEG6f);
%     EEG7a{i}=ft_freqanalysis(cfg, EEG7f);
%     EEG8a{i}=ft_freqanalysis(cfg, EEG8f);

   cfg.output        = 'original';
   EEG1b{i}=ft_freqanalysis(cfg, EEG1f);
   EEG2b{i}=ft_freqanalysis(cfg, EEG2f);
%    EEG7b{i}=ft_freqanalysis(cfg,EEG7f);
%    EEG8b{i}=ft_freqanalysis(cfg,EEG8f);
%    EEG5b{i}=ft_freqanalysis(cfg, EEG5f);
%    EEG6b{i}=ft_freqanalysis(cfg, EEG6f);

    cfg                           = [];
    cfg.parameter         = 'powspctrm';
    cfg.operation           = 'x1-x2';
    
    EEG1o{i}=ft_math(cfg,EEG1b{i}, EEG1a{i});
    EEG2o{i}=ft_math(cfg,EEG2b{i}, EEG2a{i});
%     EEG7o{i}=ft_math(cfg,EEG7b{i}, EEG7a{i});
%     EEG8o{i}=ft_math(cfg,EEG8b{i}, EEG8a{i})
%     EEG5o{i}=ft_math(cfg,EEG5b{i}, EEG5a{i});
%     EEG6o{i}=ft_math(cfg,EEG6b{i}, EEG6a{i})
    
    
    %% Coherence
    cfg                   = [];
    cfg.method     = 'mtmfft';
    cfg.taper         = 'dpss';
    cfg.output       = 'fourier';
    cfg.tapsmofrq = 2;
    
    EEG1fft         = ft_freqanalysis(cfg, EEG1f);
    EEG2fft         = ft_freqanalysis(cfg, EEG2f);
%     EEG5fft         = ft_freqanalysis(cfg, EEG5f);
%     EEG6fft         = ft_freqanalysis(cfg, EEG6f);
%     EEG7fft         = ft_freqanalysis(cfg, EEG7f);
%     EEG8fft         = ft_freqanalysis(cfg, EEG8f);

    cfg                          = [];
cfg.method                = 'coh';
coh_EEG1o {i}          = ft_connectivityanalysis(cfg, EEG1fft);
coh_EEG2o {i}          = ft_connectivityanalysis(cfg, EEG2fft);
% coh_EEG7o {i}          = ft_connectivityanalysis(cfg, EEG7fft);
% coh_EEG8o {i}          = ft_connectivityanalysis(cfg, EEG8fft);
% coh_EEG5o {i}          = ft_connectivityanalysis(cfg, EEG5fft);
% coh_EEG6o {i}          = ft_connectivityanalysis(cfg, EEG6fft);

end

%% GRAND AVERAGE
cfg=[];
cfg.keepindividual = 'yes';
 EEG1o2 = ft_freqgrandaverage(cfg,EEG1o{:});
 EEG2o2 = ft_freqgrandaverage(cfg,EEG2o{:});
%  EEG7o2 = ft_freqgrandaverage(cfg,EEG7o{:});
%  EEG8o2 = ft_freqgrandaverage(cfg,EEG8o{:});
%  EEG5o2 = ft_freqgrandaverage(cfg,EEG5o{:});
%  EEG6o2 = ft_freqgrandaverage(cfg,EEG6o{:});
%  
 EEG1f2 = ft_freqgrandaverage(cfg,EEG1a{:});
 EEG2f2 = ft_freqgrandaverage(cfg,EEG2a{:});
%  EEG7f2 = ft_freqgrandaverage(cfg,EEG7a{:});
%  EEG8f2 = ft_freqgrandaverage(cfg,EEG8a{:});
%  EEG5f2 = ft_freqgrandaverage(cfg,EEG5a{:});
%  EEG6f2 = ft_freqgrandaverage(cfg,EEG6a{:});
%  
 
% COHERENCE
 
 for i = 1:16
 coh_EEG1o2{i}= coh_EEG1o{i}; coh_EEG1o2{i}.cohspctrm = squeeze(mean(coh_EEG1o{i}.cohspctrm)); coh_EEG1o2{i}.dimord= 'chan_freq';
 coh_EEG2o2{i}= coh_EEG2o{i}; coh_EEG2o2{i}.cohspctrm = squeeze(mean(coh_EEG2o{i}.cohspctrm)); coh_EEG2o2{i}.dimord= 'chan_freq';
%  coh_EEG7o2{i}= coh_EEG7o{i}; coh_EEG7o2{i}.cohspctrm = squeeze(mean(coh_EEG7o{i}.cohspctrm)); coh_EEG7o2{i}.dimord= 'chan_freq';
%  coh_EEG8o2{i}= coh_EEG8o{i}; coh_EEG8o2{i}.cohspctrm = squeeze(mean(coh_EEG8o{i}.cohspctrm)); coh_EEG8o2{i}.dimord= 'chan_freq';
%  coh_EEG5o2{i}= coh_EEG5o{i}; coh_EEG5o2{i}.cohspctrm = squeeze(mean(coh_EEG5o{i}.cohspctrm)); coh_EEG5o2{i}.dimord= 'chan_freq';
%   coh_EEG6o2{i}= coh_EEG6o{i}; coh_EEG6o2{i}.cohspctrm = squeeze(mean(coh_EEG6o{i}.cohspctrm)); coh_EEG6o2{i}.dimord= 'chan_freq';
 end
 %% Calculate slopes
 cfg=[];
 cfg.frequency = [1 45];
%  cfg.avgoverchan = 'yes';
%  EEG7f3 = ft_selectdata(cfg,EEG7f2);
%  EEG8f3 = ft_selectdata(cfg,EEG8f2);
%  EEG5f3 = ft_selectdata(cfg,EEG5f2);
%  EEG6f3 = ft_selectdata(cfg,EEG6f2);
 EEG1f3 = ft_selectdata(cfg,EEG1f2);
 EEG2f3 = ft_selectdata(cfg,EEG2f2);

 for i = 1:16
     for ii = 1:11
% pNE = polyfit(log(EEG7f3.freq)',  squeeze(log(EEG7f3.powspctrm(i,:,:))), 1); pNE2(i,ii)=pNE(1);
% pNA= polyfit(log(EEG8f3.freq)',  squeeze(log(EEG8f3.powspctrm(i,:,:))), 1); pNA2(i,ii)=pNA(1);
% 
% pRE = polyfit(log(EEG5f3.freq)',  squeeze(log(EEG5f3.powspctrm(i,:,:))), 1); pRE2(i,ii)=pRE(1);
% pRA= polyfit(log(EEG6f3.freq)',  squeeze(log(EEG6f3.powspctrm(i,:,:))), 1); pRA2(i,ii)=pRA(1);

    pAE= polyfit(log(EEG1f3.freq)',  squeeze(log(EEG1f3.powspctrm(i,ii,:))), 1); pAE2(i,ii)=pAE(1);
    pAA = polyfit(log(EEG2f3.freq)',  squeeze(log(EEG2f3.powspctrm(i,ii,:))), 1); pAA2(i,ii)=pAA(1);
     end
 end
%  for i = 1:16
%      for ii = 1:11
% pNE = polyfit(log(EEG7f3.freq)',  squeeze(log(EEG7f3.powspctrm(i,ii,:))), 1); pNE2(i,ii)=pNE(1);
% pNA= polyfit(log(EEG8f3.freq)',  squeeze(log(EEG8f3.powspctrm(i,ii,:))), 1); pNA2(i,ii)=pNA(1);
% 
% pRE = polyfit(log(EEG5f3.freq)',  squeeze(log(EEG5f3.powspctrm(i,ii,:))), 1); pRE2(i,ii)=pRE(1);
% pRA= polyfit(log(EEG6f3.freq)',  squeeze(log(EEG6f3.powspctrm(i,ii,:))), 1); pRA2(i,ii)=pRA(1);
% 
% pAE= polyfit(log(EEG1f3.freq)',  squeeze(log(EEG1f3.powspctrm(i,ii,:))), 1); pAE2(i,ii)=pAE(1);
% pAA = polyfit(log(EEG2f3.freq)',  squeeze(log(EEG2f3.powspctrm(i,ii,:))), 1); pAA2(i,ii)=pAA(1);
%      end
%  end
 
 % ALL
sNE=pNE2(:,1);
sNA=pNA2(:,1);
% sI1=pI1(:,1);
% sI2=pI2(:,1);
sRE=pRE2(:,1);
sRA=pRA2(:,1);
 writetable(table(sNE,sNA,sRE,sRA,'VariableNames',{'NE', 'NA','RE','RA'}), 'Slope345_REM_all.txt')

 
 [p,h,stat] =  signrank(pNE2(:,3),pNA2(:,3))
 
%Performace Gain
PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);
PG_S = PG.Var1(SG); PG_W=PG.Var1(AG);


for i = 1:11
[r(i), p(i)] =corr(pAE2(:,i), PG_S, 'type','Kendall');
end



%% BETA Correlations 
cfg                = [];
cfg.keepindividual = 'yes';
cfg.frequency      = [13 30];
cfg.avgoverfreq    = 'yes';
% cfg.channel        = {'C3','C4'};
%  cfg.avgoverchan    = 'yes';
% cfg.latency        = [-1 0];
for ii = 1:16
EEG1S{ii} = ft_selectdata(cfg, EEG1o{ii});
EEG2S{ii} = ft_selectdata(cfg, EEG2o{ii});
% EEG7S{ii} = ft_selectdata(cfg, EEG7o{ii});
% EEG8S{ii} = ft_selectdata(cfg, EEG8o{ii});
% EEG5S{ii} = ft_selectdata(cfg, EEG5o{ii});
% EEG6S{ii} = ft_selectdata(cfg, EEG6o{ii});
end

A=[];B=[];C=[];D=[];E=[];F=[];
for ii = 1:16
     for iii = 1:11
% A(ii) = EEG7S{ii}.powspctrm;
% B(ii) = EEG8S{ii}.powspctrm;
% C(ii) = EEG5S{ii}.powspctrm;
% D(ii) = EEG6S{ii}.powspctrm;
E(ii,iii) = EEG1S{ii}.powspctrm(iii);
F(ii,iii) = EEG2S{ii}.powspctrm(iii);
     end
end


%C3 and C4
exp = mean(E(:,4:5),2);
adt = mean(F(:,4:5),2);

writetable(table(exp,adt, 'VariableNames', {'All_E','All_A'}), 'REMall_Beta_C.txt');



writetable(table(E',F',C',D',A',B', 'VariableNames', {'All_E','All_A','R_E','R_A','N_E','N_A'}), 'REMSleep__4s.txt');


writetable(table(E',F',C',D',A',B', 'VariableNames', {'All_E','All_A','R_E','R_A','N_E','N_A'}), 'REMSleep__4s.txt');


% Groups 
AG = [1:10, 22, 23,24,25,31,32,33];
SG = [11:21, 26,27,28,29,30];

%Performace Gain
PG = readtable('/home/b1044271/EEG_KB/Results/Behavioral/Perf_gain_all.txt', 'ReadVariableNames', 0);
PG_S = PG.Var1(SG); PG_W=PG.Var1(AG);

 [r,p] = corr(E'-F',PG_S , 'type','pearson')

 for x = 1:11
 [p(x),h,stat] = signrank (E(:,x),F(:,x))
 end
for iii = 1:11
        
[p1,h,stat] = signrank(C,D)
[p2,h,stat] = signrank(A,B)
[p3,h,stat] = signrank(E,F)

end

% 
% writetable(table(A',B','VariableNames',{'E','A'}), 'BBA_NEMs.txt');

EEG_E_SG = squeeze(mean(EEG1pow3(:,:)));
EEG_E_SGs = squeeze(std(EEG1pow3(:,:)))./sqrt(16);

EEG_A_SG = squeeze(mean(EEG2pow3(:,:)));
EEG_A_SGs = squeeze(std(EEG2pow3(:,:)))./sqrt(16);

figure;colors1 = [0 0 1; 1 0 0];
bl = boundedline(1:length(EEG_E_SG),EEG_E_SG,EEG_E_SGs, 1:length(EEG_E_SG),EEG_A_SG,EEG_A_SGs, 'cmap', colors1, 'alpha','transparency', 0.2);

[p,h,stat] = signrank(squeeze(mean(EEG1pow3,2)), squeeze(mean(EEG2pow3,2)))

EEG_E_T = squeeze(mean(EEG3pow3(:,:)));
EEG_E_Ts = squeeze(std(EEG3pow3(:,:)))./sqrt(16);

EEG_A_T = squeeze(mean(EEG4pow3(:,:)));
EEG_A_Ts = squeeze(std(EEG4pow3(:,:)))./sqrt(16);

figure;colors1 = [0 0 1; 1 0 0];
bl = boundedline(1:length(EEG_E_T(3:298)),EEG_E_T(3:298),EEG_E_Ts(3:298), 1:length(EEG_E_T(3:298)),EEG_A_T(3:298),EEG_A_Ts(3:298), 'cmap', colors1, 'alpha','transparency', 0.2);
[p,h,stat] = signrank(squeeze(nanmean(EEG3pow3,2)), squeeze(nanmean(EEG4pow3,2))) 

%% COMPARE

cfg=[];
cfg.statistic        ='ft_statfun_depsamplesT'; %between groups independant sampels statistics
cfg.method           ='montecarlo'; % Monte Carlo method for calculating the signif prob an estimate of the p-value under the permutation distribution.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
% cfg.avgoverchan      = 'yes';

 cfg.frequency        = [1 45];
cfg.channel      = {'Fz'};
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
cfg.minnbchan        = 0;
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

[st_NEM] =   ft_freqstatistics(cfg, EEG7o2 , EEG8o2);
[st_REM] =   ft_freqstatistics(cfg, EEG5o2 , EEG6o2);
[st_ALL] =   ft_freqstatistics(cfg, EEG1o2 , EEG2o2);



%% Prepapre BETA for correlations

cfg                = [];
cfg.keepindividual = 'yes';
cfg.frequency      = [13 30];
cfg.avgoverfreq    = 'yes';
cfg.channel        = {'C3','C4'};
 cfg.avgoverchan    = 'yes';
% cfg.latency        = [-1 0];
for ii = 1:16
RALL_E{ii} = ft_selectdata(cfg, EEG1o{ii});
RALL_A{ii} = ft_selectdata(cfg, EEG2o{ii});
% NEMS_E{ii} = ft_selectdata(cfg, EEG7o{ii});
% NEMS_A{ii} = ft_selectdata(cfg, EEG8o{ii});
% REMS_E{ii} = ft_selectdata(cfg, EEG5o{ii});
% REMS_A{ii} = ft_selectdata(cfg, EEG6o{ii});
end


for iii = 1:16
A = 
B = 
end

for iv = 1:16
    F3_NE(iv) = NEMS_E{iv}.powspctrm(1);
    Fz_NE(iv) = NEMS_E{iv}.powspctrm(2);
    F4_NE(iv) = NEMS_E{iv}.powspctrm(3);
    C3_NE(iv) = NEMS_E{iv}.powspctrm(4);
    C4_NE(iv) = NEMS_E{iv}.powspctrm(5);
    P3_NE(iv) = NEMS_E{iv}.powspctrm(6);
    Pz_NE(iv) = NEMS_E{iv}.powspctrm(7);
    P4_NE(iv) = NEMS_E{iv}.powspctrm(8);
    O1_NE(iv) = NEMS_E{iv}.powspctrm(9);
    O2_NE(iv) = NEMS_E{iv}.powspctrm(10);
    Cz_NE(iv) = NEMS_E{iv}.powspctrm(11);
    
    F3_NA(iv) = NEMS_A{iv}.powspctrm(1);
    Fz_NA(iv) = NEMS_A{iv}.powspctrm(2);
    F4_NA(iv) = NEMS_A{iv}.powspctrm(3);
    C3_NA(iv) = NEMS_A{iv}.powspctrm(4);
    C4_NA(iv) = NEMS_A{iv}.powspctrm(5);
    P3_NA(iv) = NEMS_A{iv}.powspctrm(6);
    Pz_NA(iv) = NEMS_A{iv}.powspctrm(7);
    P4_NA(iv) = NEMS_A{iv}.powspctrm(8);
    O1_NA(iv) = NEMS_A{iv}.powspctrm(9);
    O2_NA(iv) = NEMS_A{iv}.powspctrm(10);
    Cz_NA(iv) = NEMS_A{iv}.powspctrm(11);
    
        F3_RE(iv) = REMS_E{iv}.powspctrm(1);
    Fz_RE(iv) = REMS_E{iv}.powspctrm(2);
    F4_RE(iv) = REMS_E{iv}.powspctrm(3);
    C3_RE(iv) = REMS_E{iv}.powspctrm(4);
    C4_RE(iv) = REMS_E{iv}.powspctrm(5);
    P3_RE(iv) = REMS_E{iv}.powspctrm(6);
    Pz_RE(iv) = REMS_E{iv}.powspctrm(7);
    P4_RE(iv) = REMS_E{iv}.powspctrm(8);
    O1_RE(iv) = REMS_E{iv}.powspctrm(9);
    O2_RE(iv) = REMS_E{iv}.powspctrm(10);
    Cz_RE(iv) = REMS_E{iv}.powspctrm(11);
    
    F3_RA(iv) = REMS_A{iv}.powspctrm(1);
    Fz_RA(iv) = REMS_A{iv}.powspctrm(2);
    F4_RA(iv) = REMS_A{iv}.powspctrm(3);
    C3_RA(iv) = REMS_A{iv}.powspctrm(4);
    C4_RA(iv) = REMS_A{iv}.powspctrm(5);
    P3_RA(iv) = REMS_A{iv}.powspctrm(6);
    Pz_RA(iv) = REMS_A{iv}.powspctrm(7);
    P4_RA(iv) = REMS_A{iv}.powspctrm(8);
    O1_RA(iv) = REMS_A{iv}.powspctrm(9);
    O2_RA(iv) = REMS_A{iv}.powspctrm(10);
    Cz_RA(iv) = REMS_A{iv}.powspctrm(11);
    
    F3_AE(iv) = RALL_E{iv}.powspctrm(1);
    Fz_AE(iv) = RALL_E{iv}.powspctrm(2);
    F4_AE(iv) = RALL_E{iv}.powspctrm(3);
    C3_AE(iv) = RALL_E{iv}.powspctrm(4);
    C4_AE(iv) = RALL_E{iv}.powspctrm(5);
    P3_AE(iv) = RALL_E{iv}.powspctrm(6);
    Pz_AE(iv) = RALL_E{iv}.powspctrm(7);
    P4_AE(iv) = RALL_E{iv}.powspctrm(8);
    O1_AE(iv) = RALL_E{iv}.powspctrm(9);
    O2_AE(iv) = RALL_E{iv}.powspctrm(10);
    Cz_AE(iv) = RALL_E{iv}.powspctrm(11);
    
    F3_AA(iv) = RALL_A{iv}.powspctrm(1);
    Fz_AA(iv) = RALL_A{iv}.powspctrm(2);
    F4_AA(iv) = RALL_A{iv}.powspctrm(3);
    C3_AA(iv) = RALL_A{iv}.powspctrm(4);
    C4_AA(iv) = RALL_A{iv}.powspctrm(5);
    P3_AA(iv) = RALL_A{iv}.powspctrm(6);
    Pz_AA(iv) = RALL_A{iv}.powspctrm(7);
    P4_AA(iv) = RALL_A{iv}.powspctrm(8);
    O1_AA(iv) = RALL_A{iv}.powspctrm(9);
    O2_AA(iv) = RALL_A{iv}.powspctrm(10);
    Cz_AA(iv) = RALL_A{iv}.powspctrm(11);
end


% spindle info

% cd '/home/b1044271/EEG_KB/Sleep/Spindles_detection/spindleFeatures/spindleResults'
path_R='/home/b1044271/EEG_KB/Sleep/Spindles_detection/spindleFeatures/spindleResults/';
dataRow_r   =struct2cell(dir(fullfile(path_R , '*.xls'))');
files_R =dataRow_r(1,:);


for i = 1:length(files_R)
    
    [~,sheet_name]=xlsfinfo([path_R files_R{i}]);


    for k=1:numel(sheet_name)
      data{k}=xlsread([path_R files_R{i}],sheet_name{k});
    end

    data2 = data(~cellfun('isempty',data));
    
    %% NREM
     NAmat = data2{1};
     NAmat_A(1:16,i) = NAmat(1:2:end,15);
     NAmat_E(1:16,i)  = NAmat(2:2:end,15);

    %% N2
    N2mat = data2{2};
    N2mat_A(1:16,i)  = N2mat(1:2:end,15);
    N2mat_E(1:16,i)  = N2mat(2:2:end,15); 
    %% N3
    N3mat = data2{3};
     N3mat_A(1:16,i)  = N3mat(1:2:end,15);
     N3mat_E(1:16,i)  = N3mat(2:2:end,15);
     
end



%% correlations

 [r,p] = corr(C3_NE', N3mat_E(:,1) , 'type','Kendall')
  [r,p] = corr(C4_NE', N3mat_E(:,2) , 'type','Kendall')
  [r,p] = corr(Cz_NE', N3mat_E(:,3), 'type','Kendall') 
  [r,p] = corr(F3_NE', N3mat_E(:,4), 'type','Kendall') 
  [r,p] = corr(F4_NE', N3mat_E(:,5), 'type','Kendall')
  [r,p] = corr(Fz_NE', N3mat_E(:,6), 'type','Kendall')
  [r,p] = corr(O1_NE', N3mat_E(:,7), 'type','Kendall')
  [r,p] = corr(O2_NE', N3mat_E(:,8), 'type','Kendall')
  [r,p] = corr(P3_NE', N3mat_E(:,9), 'type','Kendall')
  [r,p] = corr(P4_NE', N3mat_E(:,10), 'type','Kendall')
  [r,p] = corr(Pz_NE', N3mat_E(:,11), 'type','Kendall')
  
   [r,p] = corr(C3_AE', PG_S, 'type','Kendall')
   [r,p] = corr(C4_AE', PG_S, 'type','Kendall')
   [r,p] = corr(Cz_RE', PG_S, 'type','Kendall')
   [r,p] = corr(F3_RE', PG_S, 'type','Kendall')
   [r,p] = corr(F4_RE', PG_S, 'type','Kendall')
   [r,p] = corr(Fz_RE', PG_S, 'type','Kendall')
   [r,p] = corr(P3_NE', PG_S, 'type','Kendall')
   [r,p] = corr(P4_NE', PG_S, 'type','Kendall')
   [r,p] = corr(Pz_NE', PG_S, 'type','Kendall')

 [r,p] = corr(C3_NE'-C3_NA', N3mat_E(:,1)-N3mat_A(:,1), 'type','Kendall')
 [r,p] = corr(C4_NE'-C4_NA', N3mat_E(:,2)-N3mat_A(:,2), 'type','Kendall')
 [r,p] = corr(Cz_NE'-Cz_NA', N3mat_E(:,3)-N3mat_A(:,3), 'type','Kendall')
 
 
 [r,p] = corr((C3_NE'./C3_RE'), N2mat_E(:,1), 'type','Kendall')
 [r,p] = corr(mean([C4_NE' ,C4_RE'],2), NAmat_E(:,2), 'type','Kendall')
 [r,p] = corr(mean([Cz_NE' ,Cz_RE'],2), NAmat_E(:,3), 'type','Kendall')
 [r,p] = corr(mean([F3_NE' ,F3_RE'],2), NAmat_E(:,4), 'type','Kendall')
 [r,p] = corr(mean([F4_NE' ,F4_RE'],2), NAmat_E(:,5), 'type','Kendall')
 [r,p] = corr(Fz_NE', N3mat_E(:,6), 'type','Kendall')
 [r,p] = corr(O1_NE', N3mat_E(:,7), 'type','Kendall')
  [r,p] = corr(O2_NE', N3mat_E(:,8), 'type','Kendall')
  [r,p] = corr(P3_NE', N3mat_E(:,9), 'type','Kendall')
  [r,p] = corr(P4_NE', N3mat_E(:,10), 'type','Kendall')
  [r,p] = corr(Pz_NE', N3mat_E(:,11), 'type','Kendall')
 
  [r,p] = corr(Cz_RE', N3mat_E(:,3), 'type','Kendall') 
  [r,p] = corr(F3_RE', N3mat_E(:,4), 'type','Kendall') 
  [r,p] = corr(F4_RE', N3mat_E(:,5), 'type','Kendall')
  [r,p] = corr(Fz_NE', N3mat_E(:,6), 'type','Kendall')
  [r,p] = corr(O1_NE', N3mat_E(:,7), 'type','Kendall')
  [r,p] = corr(O2_NE', N3mat_E(:,8), 'type','Kendall')
  [r,p] = corr(P3_NE', N3mat_E(:,9), 'type','Kendall')
  [r,p] = corr(P4_NE', N3mat_E(:,10), 'type','Kendall')
  [r,p] = corr(Pz_NE', N3mat_E(:,11), 'type','Kendall')

  
  %% Coherence permutation
  cfg=[];
  cfg.parameter        ='cohspctrm';
cfg.statistic        ='ft_statfun_depsamplesT'; %between groups independant sampels statistics
cfg.method           ='montecarlo'; % Monte Carlo method for calculating the signif prob an estimate of the p-value under the permutation distribution.
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.avgoverchan      = 'yes';
 cfg.frequency        = [1 30];
% cfg.channel      = {'Pz','P4'};
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

 load('/home/b1044271/Toolboxes/LAYOUT/ChanInfoFT.mat')
cfg.neighbours    = neighb;
[stat_c] =  ft_freqstatistics(cfg, coh_EEG5o{:},coh_EEG6o{:})
  
  [stat_c2] =  ft_freqstatistics(cfg, coh_EEG7o{:},coh_EEG8o{:})

  
  %% GET BETA COHERENCE
  cfg=[];
  cfg.frequency = [13 30];
  cfg.channel = {'C3','C4'};
  cfg.avgoverfreq = 'yes';
  for i = 1:16
  cEEG7o{i} = ft_selectdata(cfg,coh_EEG7o2{i}); cohNE(i) = mean(cEEG7o{i}.cohspctrm);
  cEEG8o{i} = ft_selectdata(cfg,coh_EEG8o2{i}); cohNA(i) = mean(cEEG8o{i}.cohspctrm);
  cEEG5o{i} = ft_selectdata(cfg,coh_EEG5o2{i}); cohRE(i) = mean(cEEG5o{i}.cohspctrm);
  cEEG6o{i} = ft_selectdata(cfg,coh_EEG6o2{i}); cohRA(i) = mean(cEEG6o{i}.cohspctrm);
  end
  
  
   writetable(table(cohNE',cohNA',cohRE',cohRA','VariableNames',{'NE', 'NA','RE','RA'}), 'Coh_Beta_REM_C_4s.txt')

  
  
