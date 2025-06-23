%%
restoredefaultpath
addpath '~/fieldtrip-20210212/'
ft_defaults

addpath('~/utils/reinstatement-analysis-master/additional_functions')
addpath('~/utils')
addpath('~/utils/functions')
addpath('~/utils/beeswarm-master')

cd ~/data2github
load TF_data.mat
load('data_hc_sametrls.mat')

subjects = [60 13 15 160 25 32 33 600 2000 36 38 8];

%%
Amylist=[2 2 2 1 1 1 2 1 1 1 2 2];
Hclist=[1 2 1 2 1 1 1 1 1 2 2 1];

freq = [23:47]
toi=[-0.5 1.5];% change toi to include baseline for R1

%% find hippo peaks during encoding and prepare the data for the ERS analysis

parfor s=1:12;
[TFpeakall_erh{s}, TFpeakall_ehith{s}, nb_er{s}] = findhcpeaksEESr1(eR_hc{s},TF_erh(s),Hclist(s), TF_ehith(s),Hclist(s),50,freq,toi);
[TFpeakall_nrh{s}, TFpeakall_nhith{s}, nb_nr{s}] = findhcpeaksEESr1(nR_hc{s},TF_nrh(s),Hclist(s), TF_nhith(s),Hclist(s),50,freq,toi);
[TFpeakall_efh{s}, TFpeakall_emissh{s}, nb_ef{s}] = findhcpeaksEESr1(eKF_hc{s},TF_efh(s),Hclist(s), TF_emissh(s),Hclist(s),50,freq,toi);
[TFpeakall_nfh{s}, TFpeakall_nmissh{s}, nb_nf{s}] = findhcpeaksEESr1(nKF_hc{s},TF_nfh(s),Hclist(s), TF_nmissh(s),Hclist(s),50,freq,toi);
end

%% compute rsa for Encoding-Encoding similarity between hippocampus and amygdala R1
nsmp    = size(TFpeakall_erh{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_ehith{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]%this is 90-150 Hz;%

parfor subj=1:12;
[allerhat{subj}] = rsa_peaks(TFpeakall_erh{subj},TFpeakall_ehith{subj},tsamps,fsh);
[allnrhat{subj}] = rsa_peaks(TFpeakall_nrh{subj},TFpeakall_nhith{subj},tsamps,fsh);
[allefhat{subj}] = rsa_peaks(TFpeakall_efh{subj},TFpeakall_emissh{subj},tsamps,fsh);
[allnfhat{subj}] = rsa_peaks(TFpeakall_nfh{subj},TFpeakall_nmissh{subj},tsamps,fsh);
end

%% changed fot hc
cd ~/data2github

load Gatemp.mat
GAs = GTFs_Multitaper_baseline_Unpl
GAs = rmfield(GAs,'powspctrm');
GAs = rmfield(GAs,'freq');

GAs.time = []
dtp= size(allerhat{1},2)-1
GAs.time = TF_erh(1).values.time(151:151+dtp);

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:151+dtp);

GAs.individual = []
Ga_erha = GAs
Ga_nrha = GAs
Ga_efha = GAs
Ga_nfha = GAs

for s=1:length(subjects)
Ga_erha.individual(s,1,:)= mean(allerhat{s},1);
Ga_nrha.individual(s,1,:)= mean(allnrhat{s},1);
Ga_efha.individual(s,1,:)= mean(allefhat{s},1);
Ga_nfha.individual(s,1,:)= mean(allnfhat{s},1);
end

%% changed for hc
URUFha = Ga_erha;
URUFha.individual = Ga_erha.individual - Ga_efha.individual
NRNFha = Ga_erha;
NRNFha.individual = Ga_nrha.individual - Ga_nfha.individual
%%
URUFhhstat = URUFha;
URUFhhstat.individual = URUFha.individual(:,1,[26:end])
URUFhhstat.time = URUFha.time(:,[26:end])
URUFhhstat.stderr = URUFha.stderr(:,[26:end])

NRNFhhstat = NRNFha;
NRNFhhstat.individual = NRNFha.individual(:,1,[26:end])
NRNFhhstat.time = NRNFha.time(:,[26:end])
NRNFhhstat.stderr = NRNFha.stderr(:,[26:end])

%% cluter stat 

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05
cfg.numrandomization = 'all';%10000;
%cfg.latency     = [0 1.5]
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(URUFha.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

int = ft_timelockstatistics(cfg,URUFhhstat, NRNFhhstat);% from 0 to 1.5

%% plot no significant results when centering on encoding hc data
%%
semeR = std(Ga_erha.individual/sqrt(ns))
semeF = std(Ga_efha.individual/sqrt(ns))
semnR = std(Ga_nrha.individual/sqrt(ns))
semnF = std(Ga_nfha.individual/sqrt(ns))

%% Supplementary Fig. 10

timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000

r=figure;set(r,'position', [0 0 1000 1500])
subplot(4,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erha.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efha.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time from hippocampus peaks (s)')
hold on
subplot(4,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrha.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfha.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time from hippocampus peaks (s)')
hold on
subplot(4,3,3)
plot(timestat,int.stat,'k','LineWidth',3)

xlabel('Time from hippocampus peaks (s)')
ylabel ('t-values')
title('emotion by memory effect');
hold on
