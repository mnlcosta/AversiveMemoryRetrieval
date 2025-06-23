
%% Fig.5
clear all
close all
clc

cd ~/data2github/
load data_ahh_ers_minus05.mat
load maxvalsimenc.mat
load('Ga_temptimeres.mat')
load 'TFtemp.mat'
load stat_ahhERS.mat

addpath ~/utils/beeswarm-master
addpath ~/utils/stdshade.m
addpath('~/utils/functions')
addpath('~/utils/reinstatement-analysis-master/additional_functions')
%%
for subj=1:12
[i,j,v]=find(maxalleR{subj} < 0)
allerahh_pos{subj}=allerahh{subj};
allerahh_pos{subj}(v,:)=[];
end

clear v
for subj=1:12
[i,j,v]=find(maxalleF{subj} < 0)
allefahh_pos{subj}=allefahh{subj};
allefahh_pos{subj}(v,:)=[];
end

clear v
for subj=1:12
[i,j,v]=find(maxallnR{subj} < 0)
allnrahh_pos{subj}=allnrahh{subj};
allnrahh_pos{subj}(v,:)=[];
end

clear v
for subj=1:12
[i,j,v]=find(maxallnF{subj} < 0)
allnfahh_pos{subj}=allnfahh{subj};
allnfahh_pos{subj}(v,:)=[];
end
%%

GAs = Ga_erah

GAs.time = []
GAs.time = TF_erh(1).values.time(151:246);

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:246);

GAs.individual = []

Ga_erahh = GAs
Ga_nrahh = GAs
Ga_efahh = GAs
Ga_nfahh = GAs

for s=1:12;
Ga_erahh.individual(s,1,:)= mean(allerahh{s},1);
Ga_nrahh.individual(s,1,:)= mean(allnrahh{s},1);
Ga_efahh.individual(s,1,:)= mean(allefahh{s},1);
Ga_nfahh.individual(s,1,:)= mean(allnfahh{s},1);
end

%%
URUFahh = Ga_erahh;
URUFahh.individual = Ga_erahh.individual - Ga_efahh.individual
NRNFahh = Ga_erahh;
NRNFahh.individual = Ga_nrahh.individual - Ga_nfahh.individual

%% cluter stat  

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05
cfg.numrandomization = 10000;
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(Ga_erahh.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

int = ft_timelockstatistics(cfg,URUFahh, NRNFahh)
%%
load stat_ahhERS.mat
%%
ns=size(Ga_erahh.individual,1);

clus =[]; i=[]; where =[]; 
    if ~isempty(int.posclusters) 
        clus = find([int.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus);
            where = find(int.posclusterslabelmat==i);
            stat_line(i,:) = [where(1) where(end)];
end

for i = 1:length(clus);
            where = find(int.negclusterslabelmat==i);
            stat_lineneg(i,:) = [where(1) where(end)];
end
%%
semeR = std(Ga_erahh.individual/sqrt(ns))
semeF = std(Ga_efahh.individual/sqrt(ns))
semnR = std(Ga_nrahh.individual/sqrt(ns))
semnF = std(Ga_nfahh.individual/sqrt(ns))

%%
meaner = squeeze(mean(Ga_erahh.individual,1))
meanef = squeeze(mean(Ga_efahh.individual,1))

timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000

r=figure;set(r,'position', [0 0 1200 1500])
subplot(4,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahh.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahh.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahh.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahh.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,3)
plot(timestat,int.stat,'k','LineWidth',3)
line([timestat(stat_line(1,1)), timestat(stat_line(1,2))],[-5, -5],'LineWidth',5, 'color', 'r','LineStyle','-'); 
ylabel ('t-values')
xlabel('Time (s)')
hold on

%% single subject plot as in Supplementary fig 16
timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000
subjects_sel = [6 13 15 16 25 32 33 6 2 36 38 8]

ns=size(allnrahh_pos,2)
r = figure;set(r,'position', [0 0 1000 600])
f=0
for v=1:12
f=f+1
subplot(4,3,f)
clear dataer dataef semer semef

dataer=(allerahh_pos{v});
dataef=(allefahh_pos{v});
semer = std(dataer/sqrt(size(dataer,1)))
semef = std(dataef/sqrt(size(dataef,1)))
shadedErrorBar(time,squeeze(mean(dataer)), semer, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(dataef)), semef, 'm', 1); hold on

ylabel('Corr Coef')
xlabel('Time (s)')
ylim([-0.5 0.5])
if f==8 | f==9 | f==12
title(['Patient Z:',num2str(subjects_sel(f))])
else
title(['Patient:',num2str(subjects_sel(f))])
end
end
%%
ns=size(allnrahh_pos,2)
r = figure;set(r,'position', [0 0 1000 600])
f=0
for v=1:12
f=f+1
subplot(4,3,f)
clear datanr datanf semnr semnf

datanr=(allnrahh_pos{v});
datanf=(allnfahh_pos{v});
semnr = std(datanr/sqrt(size(datanr,1)))
semnf = std(datanf/sqrt(size(datanf,1)))
shadedErrorBar(time,squeeze(mean(datanr)), semnr, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(datanf)), semnf, 'k', 1); hold on

ylabel('Corr Coef')
xlabel('Time (s)')
ylim([-0.5 0.5])
if f==8 | f==9 | f==12
title(['Patient Z:',num2str(subjects_sel(f))])
else
title(['Patient:',num2str(subjects_sel(f))])
end
end

