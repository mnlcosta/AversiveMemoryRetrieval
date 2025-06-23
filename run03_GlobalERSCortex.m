clear all
close all
clc

restoredefaultpath
addpath '~/fieldtrip-20210212/'
ft_defaults

cd ~/data2github
load('datapowconds_encrec_crtx.mat')
%% make subject list
subjects=[60 13 15 160 25 32 33 36 38 2000 600 8];

ts=[100:301]; %this is time until 1.5
fsh=[35:81]
win=50
%% compute rsa for each condition
for s=1:length(subjects);
all2aller = []
all2aller(:,1,1,:,:) = pow_er{s}(:,fsh,ts);
all2aller(:,2,1,:,:) = pow_ehit{s}(:,fsh,ts);

all2allekf = []
all2allekf(:,1,1,:,:) = pow_eKF{s}(:,fsh,ts);
all2allekf(:,2,1,:,:) = pow_eHitKMiss{s}(:,fsh,ts);

all2allnr = []
all2allnr(:,1,1,:,:) = pow_nr{s}(:,fsh,ts);
all2allnr(:,2,1,:,:) = pow_nhit{s}(:,fsh,ts);

all2allnkf = []
all2allnkf(:,1,1,:,:) = pow_nKF{s}(:,fsh,ts);
all2allnkf(:,2,1,:,:) = pow_nHitKMiss{s}(:,fsh,ts);

resultser{s} = rsa_m(all2aller, win, 10, [1:size(fsh,2)], 0);%all2all, win_width, mf, f, meanInTime
resultsnr{s} = rsa_m(all2allnr, win, 10, [1:size(fsh,2)], 0);
resultsekf{s} = rsa_m(all2allekf, win, 10, [1:size(fsh,2)], 0);
resultsnkf{s} = rsa_m(all2allnkf, win, 10, [1:size(fsh,2)], 0);
end
%% single subject plots
figure
for s=1:length(subjects)
subplot(3,4,s)
d2per{s} = squeeze(nanmean(resultser{s}));
imagesc(d2per{s});
title('er')
end

figure
for s=1:length(subjects)
subplot(3,4,s)
d2pnr{s} = squeeze(nanmean(resultsnr{s}));
imagesc(d2pnr{s});
title('nr')
end

figure
for s=1:length(subjects)
subplot(3,4,s)
d2pekf{s} = squeeze(nanmean(resultsekf{s}));
imagesc(d2pekf{s});
title('ekf')
end

figure
for s=1:length(subjects)
subplot(3,4,s)
d2pnkf{s} = squeeze(nanmean(resultsnkf{s}));
imagesc(d2pnkf{s});
title('nkf')
end

close all
%%
load('Gatemp.mat')

Ga_er=GTFs_Multitaper_baseline_Unpl;
Ga_er.powspctrm=[]
Ga_nr=GTFs_Multitaper_baseline_Unpl;
Ga_nr.powspctrm=[]
Ga_ekf=GTFs_Multitaper_baseline_Unpl;
Ga_ekf.powspctrm=[]
Ga_nkf=GTFs_Multitaper_baseline_Unpl;
Ga_nkf.powspctrm=[]


v=0
for s=1:length(subjects)
v=v+1
Ga_er.powspctrm(v,1,:,:)= atanh(d2per{s});
Ga_ekf.powspctrm(v,1,:,:)= atanh(d2pekf{s});
Ga_nr.powspctrm(v,1,:,:)= atanh(d2pnr{s});
Ga_nkf.powspctrm(v,1,:,:)= atanh(d2pnkf{s});
end

%% calculate diff
load('Ga_GlobalERSCrtx.mat')

URUKF = Ga_er;
URUKF.powspctrm = Ga_er.powspctrm - Ga_ekf.powspctrm

NRNKF = Ga_er;
NRNKF.powspctrm = Ga_nr.powspctrm - Ga_nkf.powspctrm

intdata = Ga_er
intdata.powspctrm = URUKF.powspctrm - NRNKF.powspctrm
%% stat int

cfg=[];
cfg.method = 'montecarlo'
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';


cfg.latency ='all';
cfg.frequency ='all';
cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 10000;

cfg.neighbours = []; 
cfg.ivar = 1;
cfg.uvar = 2;
ns=size(Ga_er.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

statERSint = ft_freqstatistics(cfg, URUKF, NRNKF);

%%
load('stat_GlobalERSCrtx.mat')

intdata = Ga_er;
intdata.powspctrm = URUKF.powspctrm - NRNKF.powspctrm
t1=-0.5103
t2=1.5009
tot=t2+abs(t1)
it=tot/16
time = [-0.5103:0.1257:1.5009]
intdata.time=[]
intdata.time=time(1,[2:end])
intdata.freq=time(1,[2:end])

%% Plot Supplementary Fig.6
cfg = [];
cfg.figure='gcf';
cfg.zlim = [-0.05 0.05]
cfg.parameter = 'powspctrm';
cfg.maskparameter='mask';
cfg.maskstyle = 'outline';

intdata.mask = logical(zeros(size(intdata.powspctrm)));
intdata.mask=(statERSint.negclusterslabelmat==1);

h = figure;set(h,'position', [0 0 1000 400]) %this is in pixel

subplot(2,3,1)
x = [-0.5 1.5];
y = [-0.5 1.5];
clims = [-0.05 0.05];
C = squeeze(mean(Ga_er.powspctrm,1))
imagesc(x,y,C,clims);
colorbar
axis('xy')
title 'eRHit'
xlabel('recognition time (s)','FontSize',10)
ylabel('encoding time (s)','FontSize',10)
hold on
plot([0,0],ylim,'k-','LineWidth',2)
hold on
y = 0;
line([-0.5,1.5],[y,y],'Color','black','LineWidth',2)
axis square

subplot(2,3,2)
x = [-0.5 1.5];
y = [-0.5 1.5];
clims = [-0.05 0.05];
C = squeeze(mean(Ga_ekf.powspctrm,1));
imagesc(x,y,C,clims);
colorbar
axis('xy')
title 'eKHitMiss'
xlabel('recognition time (s)','FontSize',10)
ylabel('encoding time (s)','FontSize',10)
hold on
plot([0,0],ylim,'k-','LineWidth',2)
hold on
y = 0;
line([-0.5,1.5],[y,y],'Color','black','LineWidth',2)
axis square

subplot(2,3,4)
x = [-0.5 1.5];
y = [-0.5 1.5];
clims = [-0.05 0.05];
C = squeeze(mean(Ga_nr.powspctrm,1))
imagesc(x,y,C,clims);
colorbar
axis('xy')
title 'nRHit'
xlabel('recognition time (s)','FontSize',10)
ylabel('encoding time (s)','FontSize',10)
hold on
plot([0,0],ylim,'k-','LineWidth',2)
hold on
y = 0;
line([-0.5,1.5],[y,y],'Color','black','LineWidth',2)
axis square

subplot(2,3,5)
x = [-0.5 1.5];
y = [-0.5 1.5];
clims = [-0.05 0.05];
C = squeeze(mean(Ga_nkf.powspctrm,1));
imagesc(x,y,C,clims);
colorbar
axis('xy')
title 'nKHitMiss'
xlabel('recognition time (s)','FontSize',10)
ylabel('encoding time (s)','FontSize',10)
hold on
plot([0,0],ylim,'k-','LineWidth',2)
hold on
y = 0;
line([-0.5,1.5],[y,y],'Color','black','LineWidth',2)
hold on
axis square

subplot(2,3,3);
cfg.xlim = [-0.5 1.5]
cfg.ylim = [-0.5 1.5]
ft_singleplotTFR(cfg,intdata);
tit=title(['ERS emotion by memory ; pval=' num2str(statERSint.posclusters(1).prob, '%0.3f')]);set(findobj(tit,'type','text'),'FontSize',10);
xlabel('recognition time (s)');ylabel('encoding time (s)');
hold on
plot([0,0],ylim,'k-','LineWidth',2)
hold on
y = 0;
line([-0.4,1.5],[y,y],'Color','black','LineWidth',2)
hold on
axis square





