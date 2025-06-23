%% This script reproduce the hippocampus time frequency results and statistic
% reported at Figure 2 and Supplementary Figure 3;

clear all
close all
clc

restoredefaultpath
addpath ' ~/fieldtrip-20210212/'
ft_defaults

cd ~/data2github/
load 'Ga_Hippocampus.mat'

subjects = [60 13 15 160 25 32 33 600 2000 36 38 8]; 

%% compute the diff for the interaction
GTFs_diff_eHitnHit = Ga_eHit
GTFs_diff_eHitnHit.powspctrm = Ga_eHit.powspctrm - Ga_nHit.powspctrm

GTFs_diff_eCrjnCrj = Ga_eCrj
GTFs_diff_eCrjnCrj.powspctrm = Ga_eCrj.powspctrm - Ga_nCrj.powspctrm

GTFs_diff_eMissnMiss = Ga_eMiss
GTFs_diff_eMissnMiss.powspctrm = Ga_eMiss.powspctrm - Ga_nMiss.powspctrm

GTFs_diff_eKMissnKMiss = Ga_eKMiss
GTFs_diff_eKMissnKMiss.powspctrm = Ga_eKMiss.powspctrm - Ga_nKMiss.powspctrm

GTFs_diff_eHitKnHitK = Ga_eHitK
GTFs_diff_eHitKnHitK.powspctrm = Ga_eHitK.powspctrm - Ga_nHitK.powspctrm
%%
cfg=[];
cfg.method = 'montecarlo'
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=35;
f2=150;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

%cfg.minnbchan      = 2;
cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 10000;

cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(GTFs_diff_eHitnHit.powspctrm,1)

cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];
Int_hitcrj = ft_freqstatistics(cfg,GTFs_diff_eHitnHit, GTFs_diff_eCrjnCrj);
Int_hitmiss = ft_freqstatistics(cfg,GTFs_diff_eHitnHit, GTFs_diff_eMissnMiss);
Int_hitKmiss = ft_freqstatistics(cfg,GTFs_diff_eHitnHit, GTFs_diff_eKMissnKMiss);

%stat for Supplementary Fig.3
Int_hithitK = ft_freqstatistics(cfg,GTFs_diff_eHitnHit, GTFs_diff_eHitKnHitK);
Int_hitmiss = ft_freqstatistics(cfg,GTFs_diff_eHitnHit, GTFs_diff_eMissnMiss);
% simple effect
eRHiteKHit = ft_freqstatistics(cfg,Ga_eHit,Ga_eHitK);
eRHiteMiss = ft_freqstatistics(cfg,Ga_eHit,Ga_eMiss);
eKHiteMiss = ft_freqstatistics(cfg,Ga_eHitK,Ga_eMiss);
%%
load statTFHippo.mat
%% plot data for Fig. 1
%% plot data for conditions
cfg = []; 
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];
cfg.ylim = [35 150];
cfg.zlim = [-0.3 0.3];
h = figure;
subplot(2,2,1) ; ft_singleplotTFR(cfg,Ga_eHit); 
tit=title('eRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,Ga_nHit);
tit=title('nRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,3); ft_singleplotTFR(cfg,Ga_eKMiss);
tit=title('eKHit&Miss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,Ga_nKMiss);
tit=title('nKHit&Miss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
%% define the cluster

h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg.zlim = [-4 4];
M = Int_hitKmiss;
M.powspctrm = Int_hitKmiss.stat.*Int_hitKmiss.mask;
M.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,M); colorbar;  tit=title('int en hit Kmiss'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16);
%% statistic plot in Fig. 2b
h = figure;
logRelative1 = Int_hitKmiss
logRelative1.mask = Int_hitKmiss.mask;
logRelative1.powspctrm = Int_hitKmiss.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);

%% mean gamma values within the significant cluster
maskt= sum(squeeze(Int_hitKmiss.mask),1);
tvec = Int_hitKmiss.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(Int_hitKmiss.mask),2);
fvec = Int_hitKmiss.freq(maskf>0);
foi=[min(fvec) max(fvec)]

t= toi
f = foi

pt1 = nearest(Ga_eHit.time,t(1));
pt2 = nearest(Ga_eHit.time,t(2));
pf1 = nearest(Ga_eHit.freq,f(1));
pf2 = nearest(Ga_eHit.freq,f(2));

clear mat 
mat(:,1) = squeeze(mean(mean(Ga_eHit.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(Ga_eKMiss.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,3) = squeeze(mean(mean(Ga_nHit.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,4) = squeeze(mean(mean(Ga_nKMiss.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
%% plot single subjects data as in Fig. 1c
addpath ~/utils/beeswarm-master
addpath ~/utils/stdshade.m

h = figure;set(h,'position', [0 0 250 188])
ns=size(mat,1)
x = [ones(ns,1) ones(ns,1)*2 ones(ns,1)*3 ones(ns,1)*4];
y = [mat(:,1) mat(:,2) mat(:,3) mat(:,4)];
beeswarm(x(:),y(:),'sort_style','up','dot_size',2,'overlay_style','sd','colormap',[1 0 0; 1 0 1;0 0 1; 0.5 0.5 0.5])
ylim([-0.6 0.8]);
ylabel('mean gamma','FontSize',10)
xticklabels({'ehit','ehitk&miss','nhit','nhitk&miss'})

%% compute the simple effect
[He,Pe,CIe,STATSe] = ttest(mat(:,1),mat(:,2));
[Hn,Pn,CIn,STATSn] = ttest(mat(:,3),mat(:,4));

 %% calculate d cohen
 der = mean(mat(:,1))
 dekf = mean(mat(:,2))
 dnr = mean(mat(:,3))
 dnkf = mean(mat(:,4))
 
 appmat_emo = [mat(:,1)-mat(:,2)];
 stdemo = std(appmat_emo)
 
 appmat_neu = [mat(:,3)-mat(:,4)];
 stdneu = std(appmat_neu)
 
 de = (der-dekf)/stdemo
 dn= (dnr-dnkf)/stdneu
%% plot data supplementary Fig. 3
cfg = []; 
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];%
cfg.ylim = [35 150];%'maxmin';
cfg.zlim = [-0.3 0.3]%'maxmin';
h = figure;%set(h,'position', [0 0 400 188])
subplot(4,2,1) ; ft_singleplotTFR(cfg,Ga_eHit); 
tit=title('eRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(4,2,2); ft_singleplotTFR(cfg,Ga_nHit);
tit=title('nRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(4,2,3) ; ft_singleplotTFR(cfg,Ga_eHitK); 
tit=title('eKHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(4,2,4); ft_singleplotTFR(cfg,Ga_nHitK);
tit=title('nKHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(4,2,5); ft_singleplotTFR(cfg,Ga_eMiss);
tit=title('eMiss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(4,2,6); ft_singleplotTFR(cfg,Ga_nMiss);
tit=title('nMiss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(4,2,7); ft_singleplotTFR(cfg,Ga_eCrj);
tit=title('eCrj');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(4,2,8); ft_singleplotTFR(cfg,Ga_nCrj);
tit=title('nCrj');set(findobj(tit,'type','text'),'FontSize',36);
%% define cluster 

h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg.zlim = [-4 4];
M = Int_hitmiss;
Int_hitmiss.mask = Int_hitmiss.posclusterslabelmat ==1 
M.powspctrm = Int_hitmiss.stat.*Int_hitmiss.mask;
M.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,M); colorbar;  tit=title('int en hit miss'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16);
%% plot statistic for Supplementary Fig. 1b
h = figure;set(h,'position', [0 0 250 188]) 

logRelative1 = Int_hitmiss
logRelative1.mask = Int_hitmiss.mask;
logRelative1.powspctrm = Int_hitmiss.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);

%% plot single subjects data as in Supplementary Fig 3c
maskt= sum(squeeze(Int_hitmiss.mask),1);
tvec = Int_hitmiss.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(Int_hitmiss.mask),2);
fvec = Int_hitmiss.freq(maskf>0);
foi=[min(fvec) max(fvec)]

t= toi
f = foi

pt1 = nearest(Ga_eHit.time,t(1));
pt2 = nearest(Ga_eHit.time,t(2));
pf1 = nearest(Ga_eHit.freq,f(1));
pf2 = nearest(Ga_eHit.freq,f(2));

clear mat 
mat(:,1) = squeeze(mean(mean(Ga_eHit.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(Ga_eMiss.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,3) = squeeze(mean(mean(Ga_nHit.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,4) = squeeze(mean(mean(Ga_nMiss.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));

addpath ~/utils/beeswarm-master
addpath ~/utils/stdshade.m

h = figure;set(h,'position', [0 0 250 188])
ns=size(mat,1)
x = [ones(ns,1) ones(ns,1)*2 ones(ns,1)*3 ones(ns,1)*4];
y = [mat(:,1) mat(:,2) mat(:,3) mat(:,4)];
beeswarm(x(:),y(:),'sort_style','up','dot_size',2,'overlay_style','sd','colormap',[1 0 0; 1 0 1;0 0 1; 0.5 0.5 0.5])
ylim([-0.6 0.8]);
ylabel('mean gamma','FontSize',10)
xticklabels({'ehit','emiss','nhit','nmiss'})

%% compute the simple effect
[He,Pe,CIe,STATSe] = ttest(mat(:,1),mat(:,2));
[Hn,Pn,CIn,STATSn] = ttest(mat(:,3),mat(:,4));

 %% calculate d cohen
 der = mean(mat(:,1))
 dekf = mean(mat(:,2))
 dnr = mean(mat(:,3))
 dnkf = mean(mat(:,4))
 
 appmat_emo = [mat(:,1)-mat(:,2)];
 stdemo = std(appmat_emo)
 
 appmat_neu = [mat(:,3)-mat(:,4)];
 stdneu = std(appmat_neu)
 
 de = (der-dekf)/stdemo
 dn= (dnr-dnkf)/stdneu
%% define the cluster for Supplementary Fig. 3d
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg.zlim = [-4 4];
M = Int_hitcrj;
Int_hitcrj.mask = Int_hitcrj.posclusterslabelmat ==1;
M.powspctrm = Int_hitcrj.stat.*Int_hitcrj.mask;
M.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,M); colorbar;  tit=title('Int hitcrj'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16);

%% plot stat 

h = figure;
logRelative1 = Int_hitcrj
logRelative1.mask = Int_hitcrj.mask;
logRelative1.powspctrm = Int_hitcrj.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);

%% plot single subjects data as in Supplementary Fig. 3d
clear mat
maskt= sum(squeeze(Int_hitcrj.mask),1);
tvec = Int_hitcrj.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(Int_hitcrj.mask),2);
fvec = Int_hitcrj.freq(maskf>0);
foi=[min(fvec) max(fvec)]

t= toi
f = foi

pt1 = nearest(Ga_eHit.time,t(1));
pt2 = nearest(Ga_eHit.time,t(2));
pf1 = nearest(Ga_eHit.freq,f(1));
pf2 = nearest(Ga_eHit.freq,f(2));

clear mat 
mat(:,1) = squeeze(mean(mean(Ga_eHit.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(Ga_eCrj.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,3) = squeeze(mean(mean(Ga_nHit.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,4) = squeeze(mean(mean(Ga_nCrj.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));

addpath ~/utils/beeswarm-master
addpath ~/utils/stdshade.m

h = figure;set(h,'position', [0 0 250 188])
ns=size(mat,1)
x = [ones(ns,1) ones(ns,1)*2 ones(ns,1)*3 ones(ns,1)*4];
y = [mat(:,1) mat(:,2) mat(:,3) mat(:,4)];
beeswarm(x(:),y(:),'sort_style','up','dot_size',2,'overlay_style','sd','colormap',[1 0 0; 1 0 1;0 0 1; 0.5 0.5 0.5])
ylim([-0.6 0.8]);
ylabel('mean gamma','FontSize',10)
xticklabels({'ehit','ecrj','nhit','ncrj'})
%%
[He,Pe,CIe,STATSe] = ttest(mat(:,1),mat(:,2));
[Hn,Pn,CIn,STATSn] = ttest(mat(:,3),mat(:,4));

 %% calculate d cohen
 der = mean(mat(:,1))
 dekf = mean(mat(:,2))
 dnr = mean(mat(:,3))
 dnkf = mean(mat(:,4))
 
 appmat_emo = [mat(:,1)-mat(:,2)];
 stdemo = std(appmat_emo)
 
 appmat_neu = [mat(:,3)-mat(:,4)];
 stdneu = std(appmat_neu)
 
 de = (der-dekf)/stdemo
 dn= (dnr-dnkf)/stdneu
%% test Low Frequency effects
%% plot low freq data as in Supplementary Fig.4
cfg = []; 
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];
cfg.ylim = [0 34];
cfg.zlim = [-1 1];
h = figure;
subplot(2,2,1) ; ft_singleplotTFR(cfg,Ga_eHit); 
tit=title('eRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,Ga_nHit);
tit=title('nRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,3); ft_singleplotTFR(cfg,Ga_eKMiss);
tit=title('eKHit&Miss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,Ga_nKMiss);
tit=title('nKHit&Miss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
%%
cfg=[];
cfg.method = 'montecarlo'
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=0;
f2=34;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

%cfg.minnbchan      = 2;
cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;%0.025;%
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 10000;

cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(GTFs_diff_eHitnHit.powspctrm,1)

cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];
Int_hitKmisslow = ft_freqstatistics(cfg,GTFs_diff_eHitnHit, GTFs_diff_eKMissnKMiss);

%%
ns=size(Ga_Unpl.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

Emolow = ft_freqstatistics(cfg,Ga_Unpl, Ga_Neu);
MemHitKMisslow= ft_freqstatistics(cfg,Ga_Hit, Ga_KMiss);

load statTFHippoLF % no significant clusters were found.
