clear all
close all
clc

restoredefaultpath
addpath '~/fieldtrip-20210212/'
ft_defaults

addpath('~/utils/reinstatement-analysis-master/additional_functions')
addpath('~/utils')
addpath('~/utils/functions')
addpath('~/utils/beeswarm-master')
%%

cd ~/data2github
load TF_data.mat
load('data_amyhc_encrec.mat')

subjects = [60 13 15 160 25 32 33 600 2000 36 38 8];

%%
Amylist=[2 2 2 1 1 1 2 1 1 1 2 2];
Hclist=[1 2 1 2 1 1 1 1 1 2 2 1];

freq = [23:47]
toi=[-0.5 1.5];% toi includes baseline

%% find amygdala peaks during encoding and prepare the data for the EES and ERS analysis
parfor s=1:12;
[TFpeakall_era{s}, TFpeakall_erh{s}, TFpeakall_ehith{s}, nb_er{s}] = findamypeaksEESERSr1(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),TF_ehith(s),Hclist(s),50,freq,toi);
[TFpeakall_nra{s}, TFpeakall_nrh{s}, TFpeakall_nhith{s}, nb_nr{s}] = findamypeaksEESERSr1(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),TF_nhith(s),Hclist(s),50,freq,toi);
[TFpeakall_efa{s}, TFpeakall_efh{s}, TFpeakall_emissh{s}, nb_ef{s}] = findamypeaksEESERSr1(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),TF_emissh(s),Hclist(s),50,freq,toi);
[TFpeakall_nfa{s}, TFpeakall_nfh{s}, TFpeakall_nmissh{s}, nb_nf{s}] = findamypeaksEESERSr1(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),TF_nmissh(s),Hclist(s),50,freq,toi);
end
%% control analysis 300 ipd Supplementary Fig. 13 & 14
parfor s=1:12;
[TFpeakall_erac1{s}, TFpeakall_erhc1{s}, TFpeakall_ehithc1{s}, nb_er{s}] = findamypeaksEESERSr1(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),TF_ehith(s),Hclist(s),150,freq,toi);
[TFpeakall_nrac1{s}, TFpeakall_nrhc1{s}, TFpeakall_nhithc1{s}, nb_nr{s}] = findamypeaksEESERSr1(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),TF_nhith(s),Hclist(s),150,freq,toi);
[TFpeakall_efac1{s}, TFpeakall_efhc1{s}, TFpeakall_emisshc1{s}, nb_ef{s}] = findamypeaksEESERSr1(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),TF_emissh(s),Hclist(s),150,freq,toi);
[TFpeakall_nfac1{s}, TFpeakall_nfhc1{s}, TFpeakall_nmisshc1{s}, nb_nf{s}] = findamypeaksEESERSr1(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),TF_nmissh(s),Hclist(s),150,freq,toi);
end
%% control analysis random peak selection Supplementary Fig. 13 & 14
parfor s=1:12;
[TFpeakall_erarnd{s}, TFpeakall_erhrnd{s}, TFpeakall_ehithrnd{s}] = findamypeaksEESERS_rndr1(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),TF_ehith(s),Hclist(s),50,freq,toi);
[TFpeakall_nrarnd{s}, TFpeakall_nrhrnd{s}, TFpeakall_nhithrnd{s}] = findamypeaksEESERS_rndr1(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),TF_nhith(s),Hclist(s),50,freq,toi);
[TFpeakall_efarnd{s}, TFpeakall_efhrnd{s}, TFpeakall_emisshrnd{s}] = findamypeaksEESERS_rndr1(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),TF_emissh(s),Hclist(s),50,freq,toi);
[TFpeakall_nfarnd{s}, TFpeakall_nfhrnd{s}, TFpeakall_nmisshrnd{s}] = findamypeaksEESERS_rndr1(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),TF_nmissh(s),Hclist(s),50,freq,toi);
end
%% for R1 test early and late peaks Supplementary Fig. 15
parfor s=1:12;
[TFpeakall_era_e{s}, TFpeakall_era_l{s}, TFpeakall_erh_e{s}, TFpeakall_erh_l{s}, nb_er{s}] = findamypeaksEES_earlylater1(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),Hclist(s),50,freq,toi);
[TFpeakall_nra_e{s}, TFpeakall_nra_l{s}, TFpeakall_nrh_e{s}, TFpeakall_nrh_l{s}, nb_nr{s}] = findamypeaksEES_earlylater1(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),Hclist(s),50,freq,toi);
[TFpeakall_efa_e{s}, TFpeakall_efa_l{s}, TFpeakall_efh_e{s}, TFpeakall_efh_l{s}, nb_ef{s}] = findamypeaksEES_earlylater1(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),Hclist(s),50,freq,toi);
[TFpeakall_nfa_e{s}, TFpeakall_nfa_l{s}, TFpeakall_nfh_e{s}, TFpeakall_nfh_l{s}, nb_nf{s}] = findamypeaksEES_earlylater1(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),Hclist(s),50,freq,toi);
end
%% compute rsa for Encoding-Encoding similarity between amygdala and hippocampus
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_erh{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]%this is 90-150 Hz;%

parfor subj=1:12;
[alleraht{subj}] = rsa_peaks(TFpeakall_era{subj},TFpeakall_erh{subj},tsamps,fsh);
[allnraht{subj}] = rsa_peaks(TFpeakall_nra{subj},TFpeakall_nrh{subj},tsamps,fsh);
[allefaht{subj}] = rsa_peaks(TFpeakall_efa{subj},TFpeakall_efh{subj},tsamps,fsh);
[allnfaht{subj}] = rsa_peaks(TFpeakall_nfa{subj},TFpeakall_nfh{subj},tsamps,fsh);
end

parfor subj=1:12;
[alleraht_c1{subj}] = rsa_peaks(TFpeakall_erac1{subj},TFpeakall_erhc1{subj},tsamps,fsh);
[allnraht_c1{subj}] = rsa_peaks(TFpeakall_nrac1{subj},TFpeakall_nrhc1{subj},tsamps,fsh);
[allefaht_c1{subj}] = rsa_peaks(TFpeakall_efac1{subj},TFpeakall_efhc1{subj},tsamps,fsh);
[allnfaht_c1{subj}] = rsa_peaks(TFpeakall_nfac1{subj},TFpeakall_nfhc1{subj},tsamps,fsh);
end

for subj=1:12;
[alleraht_rnd{subj}] = rsa_peaks(TFpeakall_erarnd{subj},TFpeakall_erhrnd{subj},tsamps,fsh);
[allnraht_rnd{subj}] = rsa_peaks(TFpeakall_nrarnd{subj},TFpeakall_nrhrnd{subj},tsamps,fsh);
[allefaht_rnd{subj}] = rsa_peaks(TFpeakall_efarnd{subj},TFpeakall_efhrnd{subj},tsamps,fsh);
[allnfaht_rnd{subj}] = rsa_peaks(TFpeakall_nfarnd{subj},TFpeakall_nfhrnd{subj},tsamps,fsh);
end
%% compute rsa early peaks EES amy-hc
nsmp    = size(TFpeakall_era_e{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_erh_e{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]%this is 90-150 Hz;%

parfor subj=1:12;
[alleraht_e{subj}] = rsa_peaks(TFpeakall_era_e{subj},TFpeakall_erh_e{subj},tsamps,fsh);
[allnraht_e{subj}] = rsa_peaks(TFpeakall_nra_e{subj},TFpeakall_nrh_e{subj},tsamps,fsh);
[allefaht_e{subj}] = rsa_peaks(TFpeakall_efa_e{subj},TFpeakall_efh_e{subj},tsamps,fsh);
[allnfaht_e{subj}] = rsa_peaks(TFpeakall_nfa_e{subj},TFpeakall_nfh_e{subj},tsamps,fsh);
end

%% compute rsa late peaks
nsmp    = size(TFpeakall_era_l{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_erh_l{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]%this is 90-150 Hz;%

parfor subj=1:12;
[alleraht_l{subj}] = rsa_peaks(TFpeakall_era_l{subj},TFpeakall_erh_l{subj},tsamps,fsh);
[allnraht_l{subj}] = rsa_peaks(TFpeakall_nra_l{subj},TFpeakall_nrh_l{subj},tsamps,fsh);
[allefaht_l{subj}] = rsa_peaks(TFpeakall_efa_l{subj},TFpeakall_efh_l{subj},tsamps,fsh);
[allnfaht_l{subj}] = rsa_peaks(TFpeakall_nfa_l{subj},TFpeakall_nfh_l{subj},tsamps,fsh);
end
%% compute rsa for Encoding-Retrieval similarity between amygdala and hippocampus
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_ehith{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]

parfor subj=1:12;
[alleraht_ers{subj}] = rsa_peaks(TFpeakall_era{subj},TFpeakall_ehith{subj},tsamps,fsh);
[allnraht_ers{subj}] = rsa_peaks(TFpeakall_nra{subj},TFpeakall_nhith{subj},tsamps,fsh);
[allefaht_ers{subj}] = rsa_peaks(TFpeakall_efa{subj},TFpeakall_emissh{subj},tsamps,fsh);
[allnfaht_ers{subj}] = rsa_peaks(TFpeakall_nfa{subj},TFpeakall_nmissh{subj},tsamps,fsh);
end

parfor subj=1:12;
[alleraht_c1_ers{subj}] = rsa_peaks(TFpeakall_erac1{subj},TFpeakall_ehithc1{subj},tsamps,fsh);
[allnraht_c1_ers{subj}] = rsa_peaks(TFpeakall_nrac1{subj},TFpeakall_nhithc1{subj},tsamps,fsh);
[allefaht_c1_ers{subj}] = rsa_peaks(TFpeakall_efac1{subj},TFpeakall_emisshc1{subj},tsamps,fsh);
[allnfaht_c1_ers{subj}] = rsa_peaks(TFpeakall_nfac1{subj},TFpeakall_nmisshc1{subj},tsamps,fsh);
end

for subj=1:12;
[alleraht_rnd_ers{subj}] = rsa_peaks(TFpeakall_erarnd{subj},TFpeakall_ehithrnd{subj},tsamps,fsh);
[allnraht_rnd_ers{subj}] = rsa_peaks(TFpeakall_nrarnd{subj},TFpeakall_nhithrnd{subj},tsamps,fsh);
[allefaht_rnd_ers{subj}] = rsa_peaks(TFpeakall_efarnd{subj},TFpeakall_emisshrnd{subj},tsamps,fsh);
[allnfaht_rnd_ers{subj}] = rsa_peaks(TFpeakall_nfarnd{subj},TFpeakall_nmisshrnd{subj},tsamps,fsh);
end

%% To reproduce the figures run from here. 
subjects = [60 13 15 160 25 32 33 600 2000 36 38 8];

load Gatemp.mat 
load('dataEES_minus05.mat')
load('dataERS_minus05.mat')
load('TF_data.mat')

GAs = GTFs_Multitaper_baseline_Unpl
GAs = rmfield(GAs,'powspctrm');
GAs = rmfield(GAs,'freq');

GAs.time = []
dtp= size(alleraht{1},2)-1
GAs.time = TF_erh(1).values.time(151:151+dtp);

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:151+dtp);


GAs.individual = []
Ga_erah = GAs
Ga_nrah = GAs
Ga_efah = GAs
Ga_nfah = GAs


for s=1:length(subjects)
Ga_erah.individual(s,1,:)= mean(alleraht{s},1);
Ga_nrah.individual(s,1,:)= mean(allnraht{s},1);
Ga_efah.individual(s,1,:)= mean(allefaht{s},1);
Ga_nfah.individual(s,1,:)= mean(allnfaht{s},1);
end

%% early late Supplementary fig. 15
Ga_erah_e = GAs
Ga_nrah_e = GAs
Ga_efah_e = GAs
Ga_nfah_e = GAs


for s=1:length(subjects)
Ga_erah_e.individual(s,1,:)= mean(alleraht_e{s},1);
Ga_nrah_e.individual(s,1,:)= mean(allnraht_e{s},1);
Ga_efah_e.individual(s,1,:)= mean(allefaht_e{s},1);
Ga_nfah_e.individual(s,1,:)= mean(allnfaht_e{s},1);
end

Ga_erah_l = GAs
Ga_nrah_l = GAs
Ga_efah_l = GAs
Ga_nfah_l = GAs


for s=1:length(subjects)
Ga_erah_l.individual(s,1,:)= mean(alleraht_l{s},1);
Ga_nrah_l.individual(s,1,:)= mean(allnraht_l{s},1);
Ga_efah_l.individual(s,1,:)= mean(allefaht_l{s},1);
Ga_nfah_l.individual(s,1,:)= mean(allnfaht_l{s},1);
end

%% control analysis 
Ga_erahc1 = GAs
Ga_nrahc1 = GAs
Ga_efahc1 = GAs
Ga_nfahc1 = GAs

for s=1:length(subjects)
Ga_erahc1.individual(s,1,:)= mean(alleraht_c1{s},1);
Ga_nrahc1.individual(s,1,:)= mean(allnraht_c1{s},1);
Ga_efahc1.individual(s,1,:)= mean(allefaht_c1{s},1);
Ga_nfahc1.individual(s,1,:)= mean(allnfaht_c1{s},1);
end

%% control analysis rnd sel
Ga_erah_rnd = GAs
Ga_nrah_rnd = GAs
Ga_efah_rnd = GAs
Ga_nfah_rnd = GAs

for s=1:length(subjects)
Ga_erah_rnd.individual(s,1,:)= mean(alleraht_rnd{s},1);
Ga_nrah_rnd.individual(s,1,:)= mean(allnraht_rnd{s},1);
Ga_efah_rnd.individual(s,1,:)= mean(allefaht_rnd{s},1);
Ga_nfah_rnd.individual(s,1,:)= mean(allnfaht_rnd{s},1);
end

%% compute the diff for the interaction
URUFah = Ga_erah;
URUFah.individual = Ga_erah.individual - Ga_efah.individual
NRNFah = Ga_erah;
NRNFah.individual = Ga_nrah.individual - Ga_nfah.individual

URUFahc1 = Ga_erahc1;
URUFahc1.individual = Ga_erahc1.individual - Ga_efahc1.individual
NRNFahc1 = Ga_erahc1;
NRNFahc1.individual = Ga_nrahc1.individual - Ga_nfahc1.individual

URUFah_rnd = Ga_erah_rnd;
URUFah_rnd.individual = Ga_erah_rnd.individual - Ga_efah_rnd.individual
NRNFah_rnd = Ga_erah_rnd;
NRNFah_rnd.individual = Ga_nrah_rnd.individual - Ga_nfah_rnd.individual

URUFah_e = Ga_erah_e;
URUFah_e.individual = Ga_erah_e.individual - Ga_efah_e.individual
NRNFah_e = Ga_erah_e;
NRNFah_e.individual = Ga_nrah_e.individual - Ga_nfah_e.individual

URUFah_l = Ga_erah_l;
URUFah_l.individual = Ga_erah_l.individual - Ga_efah_l.individual
NRNFah_l = Ga_erah_l;
NRNFah_l.individual = Ga_nrah_l.individual - Ga_nfah_l.individual
%% select time from 0 to 1.5
URUFah_e_stat = Ga_erah_e;
URUFah_e_stat.individual = URUFah_e.individual(:,1,[26:end])
URUFah_e_stat.time = URUFah_e.time(1,[26:end])
URUFah_e_stat.stderr = URUFah_e.stderr(1,[26:end])

NRNFah_e_stat = Ga_erah_e;
NRNFah_e_stat.individual = NRNFah_e.individual(:,1,[26:end])
NRNFah_e_stat.time = NRNFah_e.time(1,[26:end])
NRNFah_e_stat.stderr = NRNFah_e.stderr(1,[26:end])


URUFah_l_stat = Ga_erah_l;
URUFah_l_stat.individual = URUFah_l.individual(:,1,[26:end])
URUFah_l_stat.time = URUFah_l.time(1,[26:end])
URUFah_l_stat.stderr = URUFah_l.stderr(1,[26:end])

NRNFah_l_stat = Ga_erah_l;
NRNFah_l_stat.individual = NRNFah_l.individual(:,1,[26:end])
NRNFah_l_stat.time = NRNFah_l.time(1,[26:end])
NRNFah_l_stat.stderr = NRNFah_l.stderr(1,[26:end])
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
%cfg.latency     = [0 1.5]
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(URUFah.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

URUFstat = ft_timelockstatistics(cfg,Ga_erah,Ga_efah)
NRNFstat = ft_timelockstatistics(cfg,Ga_nrah,Ga_nfah)

int_e = ft_timelockstatistics(cfg,URUFah_e_stat, NRNFah_e_stat)
int_l = ft_timelockstatistics(cfg,URUFah_l_stat, NRNFah_l_stat)

int = ft_timelockstatistics(cfg,URUFah, NRNFah)
intc1= ft_timelockstatistics(cfg,URUFahc1, NRNFahc1)
intrnd= ft_timelockstatistics(cfg,URUFah_rnd, NRNFah_rnd)

%% ERS results
GAs.time = []
GAs.time = TF_erh(1).values.time(151:151+dtp)

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:151+dtp)

GAs.individual = []
Ga_erahers = GAs
Ga_nrahers = GAs
Ga_efahers = GAs
Ga_nfahers = GAs

for s=1:length(subjects)
Ga_erahers.individual(s,1,:)= mean(alleraht_ers{s},1);
Ga_nrahers.individual(s,1,:)= mean(allnraht_ers{s},1);
Ga_efahers.individual(s,1,:)= mean(allefaht_ers{s},1);
Ga_nfahers.individual(s,1,:)= mean(allnfaht_ers{s},1);
end


Ga_erahc1ers = GAs
Ga_nrahc1ers = GAs
Ga_efahc1ers = GAs
Ga_nfahc1ers = GAs

for s=1:length(subjects)
Ga_erahc1ers.individual(s,1,:)= mean(alleraht_c1_ers{s},1);
Ga_nrahc1ers.individual(s,1,:)= mean(allnraht_c1_ers{s},1);
Ga_efahc1ers.individual(s,1,:)= mean(allefaht_c1_ers{s},1);
Ga_nfahc1ers.individual(s,1,:)= mean(allnfaht_c1_ers{s},1);
end


Ga_erah_rnders = GAs
Ga_nrah_rnders = GAs
Ga_efah_rnders = GAs
Ga_nfah_rnders = GAs

for s=1:length(subjects)
Ga_erah_rnders.individual(s,1,:)= mean(alleraht_rnd_ers{s},1);
Ga_nrah_rnders.individual(s,1,:)= mean(allnraht_rnd_ers{s},1);
Ga_efah_rnders.individual(s,1,:)= mean(allefaht_rnd_ers{s},1);
Ga_nfah_rnders.individual(s,1,:)= mean(allnfaht_rnd_ers{s},1);
end

%% 
URUFahers = Ga_erahers;
URUFahers.individual = Ga_erahers.individual - Ga_efahers.individual
NRNFahers = Ga_erahers;
NRNFahers.individual = Ga_nrahers.individual - Ga_nfahers.individual


URUFahc1ers = Ga_erahc1ers;
URUFahc1ers.individual = Ga_erahc1ers.individual - Ga_efahc1ers.individual
NRNFahc1ers = Ga_erahc1ers;
NRNFahc1ers.individual = Ga_nrahc1ers.individual - Ga_nfahc1ers.individual


URUFah_rnders = Ga_erah_rnders;
URUFah_rnders.individual = Ga_erah_rnders.individual - Ga_efah_rnders.individual
NRNFah_rnders = Ga_erah_rnders;
NRNFah_rnders.individual = Ga_nrah_rnders.individual - Ga_nfah_rnders.individual

ns=size(URUFahers.individual,1)

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
%cfg.latency     = [0 1.5]
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(URUFahers.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

URUFersstat = ft_timelockstatistics(cfg,Ga_erahers,Ga_efahers)
NRNFersstat = ft_timelockstatistics(cfg,Ga_nrahers,Ga_nfahers)

int_ers = ft_timelockstatistics(cfg,URUFahers, NRNFahers)
intc1_ers= ft_timelockstatistics(cfg,URUFahc1ers, NRNFahc1ers)
intrnd_ers= ft_timelockstatistics(cfg,URUFah_rnders, NRNFah_rnders)
%%
timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000

ns=size(Ga_erah.individual,1);

clus =[]; i=[]; where =[]; 
    if ~isempty(int.posclusters) 
        clus = find([int.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus);
            where = find(int.posclusterslabelmat==i);
            stat_line(i,:) = [where(1) where(end)];
end

clus =[]; i=[]; where =[]; 
    if ~isempty(int_ers.posclusters) 
        clus = find([int_ers.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus);
            where = find(int_ers.posclusterslabelmat==i);
            stat_lineers(i,:) = [where(1) where(end)];
end

%% to plot results in Fig.4 and Supplementary Figs


restoredefaultpath
addpath '~/fieldtrip-20210212/'
ft_defaults

addpath('~/utils/reinstatement-analysis-master/additional_functions')
addpath('~/utils/functions')
addpath('~/utils/beeswarm-master')
addpath('~/utils')

%%
cd ~/data2github
load  statEES.mat;
load('statEES_earlylate.mat')
load  statERS.mat
%%

clus =[]; i=[]; where =[]; 
    if ~isempty(int_ers.posclusters) 
        clus = find([int_ers.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus);
            where = find(int_ers.posclusterslabelmat==i);
            stat_lineers(i,:) = [where(1) where(end)];
end

clusc1 =[]; i=[]; wherec1 =[]; 
    if ~isempty(intc1_ers.posclusters) 
        clusc1 = find([intc1_ers.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clusc1);
            wherec1 = find(intc1_ers.posclusterslabelmat==i);
            stat_linec1_ers(i,:) = [wherec1(1) wherec1(end)];
end

clus_rnd =[]; i=[]; where_rnd =[]; 
    if ~isempty(intrnd_ers.posclusters) 
        clus_rnd = find([intrnd_ers.posclusters.prob]'<0.08);% result is marked with a different color
    end 
    
for i = 1:length(clus_rnd);
            where_rnd = find(intrnd_ers.posclusterslabelmat==i);
            stat_line_rnd_ers(i,:) = [where_rnd(1) where_rnd(end)];
end
%%
semeRers = std(Ga_erahers.individual/sqrt(ns))
semeFers = std(Ga_efahers.individual/sqrt(ns))
semnRers = std(Ga_nrahers.individual/sqrt(ns))
semnFers = std(Ga_nfahers.individual/sqrt(ns))

semeRc1 = std(Ga_erahc1ers.individual/sqrt(ns))
semeFc1 = std(Ga_efahc1ers.individual/sqrt(ns))
semnRc1 = std(Ga_nrahc1ers.individual/sqrt(ns))
semnFc1 = std(Ga_nfahc1ers.individual/sqrt(ns))

semeRrnd = std(Ga_erah_rnders.individual/sqrt(ns))
semeFrnd = std(Ga_efah_rnders.individual/sqrt(ns))
semnRrnd = std(Ga_nrah_rnders.individual/sqrt(ns))
semnFrnd = std(Ga_nfah_rnders.individual/sqrt(ns))
%%
ns=size(Ga_erah.individual,1);

clus =[]; i=[]; where =[]; 
    if ~isempty(int.posclusters) 
        clus = find([int.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus);
            where = find(int.posclusterslabelmat==i);
            stat_line(i,:) = [where(1) where(end)];
end

clusc1 =[]; i=[]; wherec1 =[]; 
    if ~isempty(intc1.posclusters) 
        clusc1 = find([intc1.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clusc1);
            wherec1 = find(intc1.posclusterslabelmat==i);
            stat_linec1(i,:) = [wherec1(1) wherec1(end)];
end

clus_rnd =[]; i=[]; where_rnd =[]; 
    if ~isempty(intrnd.posclusters) 
        clus_rnd = find([intrnd.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus_rnd);
            where_rnd = find(intrnd.posclusterslabelmat==i);
            stat_line_rnd(i,:) = [where_rnd(1) where_rnd(end)];
end

semeR = std(Ga_erah.individual/sqrt(ns))
semeF = std(Ga_efah.individual/sqrt(ns))
semnR = std(Ga_nrah.individual/sqrt(ns))
semnF = std(Ga_nfah.individual/sqrt(ns))

semeRc1 = std(Ga_erahc1.individual/sqrt(ns))
semeFc1 = std(Ga_efahc1.individual/sqrt(ns))
semnRc1 = std(Ga_nrahc1.individual/sqrt(ns))
semnFc1 = std(Ga_nfahc1.individual/sqrt(ns))

semeRrnd = std(Ga_erah_rnd.individual/sqrt(ns))
semeFrnd = std(Ga_efah_rnd.individual/sqrt(ns))
semnRrnd = std(Ga_nrah_rnd.individual/sqrt(ns))
semnFrnd = std(Ga_nfah_rnd.individual/sqrt(ns))
%% plot Fig. 4
timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000

r = figure;set(r,'position', [0 0 1200 400])
subplot(2,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahers.individual,1)), semeRers, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahers.individual,1)), semeFers, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
hold on
subplot(2,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahers.individual,1)), semnRers, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahers.individual,1)), semnFers, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
hold on
subplot(2,3,3)
plot(timestat,int_ers.stat,'k','LineWidth',3)
line([timestat(stat_lineers(1,1)), timestat(stat_lineers(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); %[0.9 0.5 0.6]
hold on
line([timestat(stat_lineers(2,1)), timestat(stat_lineers(2,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('hc time from amy peaks (s)')
ylabel ('t-values')
title('emotion by memory: ERS amy hc')
hold on
subplot(2,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erah.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efah.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
hold on
subplot(2,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrah.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfah.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
hold on
subplot(2,3,6)
plot(timestat,int.stat,'k','LineWidth',3)
line([timestat(stat_line(1,1)), timestat(stat_line(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('hc time from amy peaks (s)')
ylabel ('t-values')
title('emotion by memory: EES amy hc')

%% plot single subject with trial variability
%% EES results Supplementary Fig. 11
load Gatemp.mat
GAs = GTFs_Multitaper_baseline_Unpl
GAs = rmfield(GAs,'powspctrm');
GAs = rmfield(GAs,'freq');

% eRHit, eKHit&Miss 
timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000
subjects_sel = [6 13 15 16 25 32 33 6 2 36 38 8]


r = figure;set(r,'position', [0 0 1000 600])
f=0
for v=1:12
f=f+1
subplot(4,3,f)
clear dataer dataef semer semef
dataer=(alleraht{v});
dataef=(allefaht{v});
semer = std(dataer/sqrt(size(dataer,1)))
semef = std(dataef/sqrt(size(dataef,1)))
shadedErrorBar(time,squeeze(mean(dataer)), semer, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(dataef)), semef, 'm', 1); hold on

ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
ylim([-0.25 0.25])
if f==8 | f==9 | f==12
title(['Patient Z:',num2str(subjects_sel(f))])
else
title(['Patient:',num2str(subjects_sel(f))])
end
end
% nRHit, nKHit&Miss 
ns=size(subjects,2)
r = figure;set(r,'position', [0 0 1000 600])
f=0
for v=1:12
f=f+1
subplot(4,3,f)
clear datanr datanf semnr semnf

datanr=(allnraht{v});
datanf=(allnfaht{v});
semnr = std(datanr/sqrt(size(datanr,1)))
semnf = std(datanf/sqrt(size(datanf,1)))
shadedErrorBar(time,squeeze(mean(datanr)), semnr, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(datanf)), semnf, 'k', 1); hold on

ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
ylim([-0.25 0.25])
if f==8 | f==9 | f==12
title(['Patient Z:',num2str(subjects_sel(f))])
else
title(['Patient:',num2str(subjects_sel(f))])
end
end
%% Supplementary Fig. 9 ERS single subjects
% eRHit, eKHit&Miss 
ns=size(subjects,2)
r = figure;set(r,'position', [0 0 1000 600])
f=0
for v=1:12
f=f+1
subplot(4,3,f)
clear dataer dataef semer semef

dataer=(alleraht_ers{v});
dataef=(allefaht_ers{v});
semer = std(dataer/sqrt(size(dataer,1)))
semef = std(dataef/sqrt(size(dataef,1)))
shadedErrorBar(time,squeeze(mean(dataer)), semer, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(dataef)), semef, 'm', 1); hold on

ylabel('Corr Coef')
xlabel('Time (s)')
ylim([-0.25 0.3])
if f==8 | f==9 | f==12
title(['Patient Z:',num2str(subjects_sel(f))])
else
title(['Patient:',num2str(subjects_sel(f))])
end
end
% nRHit, nKHit&Miss 
ns=size(subjects,2)
r = figure;set(r,'position', [0 0 1000 600])
f=0
for v=1:12
f=f+1
subplot(4,3,f)
clear datanr datanf semnr semnf

datanr=(allnraht_ers{v});
datanf=(allnfaht_ers{v});
semnr = std(datanr/sqrt(size(datanr,1)))
semnf = std(datanf/sqrt(size(datanf,1)))
shadedErrorBar(time,squeeze(mean(datanr)), semnr, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(datanf)), semnf, 'k', 1); hold on

ylabel('Corr Coef')
xlabel('Time (s)')
ylim([-0.25 0.25])
if f==8 | f==9 | f==12
title(['Patient Z:',num2str(subjects_sel(f))])
else
title(['Patient:',num2str(subjects_sel(f))])
end
end

%% plot Supplementary Fig.13
timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000

r=figure;set(r,'position', [0 0 1200 1500])
subplot(4,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahc1ers.individual,1)), semeRc1, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahc1ers.individual,1)), semeFc1, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahc1ers.individual,1)), semnRc1, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahc1ers.individual,1)), semnFc1, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,3)
plot(timestat,intc1_ers.stat,'k','LineWidth',3)
line([timestat(stat_linec1_ers(1,1)), timestat(stat_linec1_ers(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
hold on
line([timestat(stat_linec1_ers(2,1)), timestat(stat_linec1_ers(2,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('Time (s)')
ylabel ('t-values')
title(['Emotion by memory int; inp=300 p=' num2str(intc1_ers.posclusters(1).prob,'%0.3f')])
hold on
subplot(4,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erah_rnders.individual,1)), semeRrnd, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efah_rnders.individual,1)), semeFrnd, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrah_rnders.individual,1)), semnRrnd, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfah_rnders.individual,1)), semnFrnd, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,6)
plot(timestat,intrnd_ers.stat,'k','LineWidth',3)
line([timestat(stat_line_rnd_ers(1,1)), timestat(stat_line_rnd_ers(1,2))],[-2, -2],'LineWidth',5, 'color', [0.9 0.5 0.6],'LineStyle','-'); 
xlabel('Time (s)')
ylabel ('t-values')
title(['Emotion by memory int; rnd peak p=' num2str(intrnd_ers.posclusters(1).prob, '%0.3f')])

%% Supplementary Fig. 14
%% plot
timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000

r=figure;set(r,'position', [0 0 1200 1500])
subplot(4,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahc1.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahc1.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time from amygdala peaks (s)')
hold on
subplot(4,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahc1.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahc1.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time from amygdala peaks (s)')
hold on
subplot(4,3,3)
plot(timestat,intc1.stat,'k','LineWidth',3)
line([timestat(stat_linec1(1,1)), timestat(stat_linec1(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('Time from amygdala peaks (s)')
ylabel ('t-values')
title('emotion by memory effect');
hold on
subplot(4,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erah_rnd.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efah_rnd.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time from amygdala peaks (s)')
hold on
subplot(4,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrah_rnd.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfah_rnd.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time from amygdala peaks (s)')
hold on
subplot(4,3,6)
plot(timestat,intrnd.stat,'k','LineWidth',3)
line([timestat(stat_line_rnd(1,1)), timestat(stat_line_rnd(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('Time from amygdala peaks (s)')
ylabel ('t-values')
title('emotion by memory effect');
hold on
%% plot results early late; Supplementary Fig. 15
ns=size(Ga_erah_e.individual,1);

clus =[]; i=[]; where =[]; 
    if ~isempty(int_l.posclusters) 
        clus = find([int_l.posclusters.prob]'<0.05);
    end 
    
for i = 1:length(clus);
            where = find(int_l.posclusterslabelmat==i);
            stat_line_l(i,:) = [where(1) where(end)];
end

semeR_l = std(Ga_erah_l.individual/sqrt(ns))
semeF_l = std(Ga_efah_l.individual/sqrt(ns))
semnR_l = std(Ga_nrah_l.individual/sqrt(ns))
semnF_l = std(Ga_nfah_l.individual/sqrt(ns))

semeR_e = std(Ga_erah_e.individual/sqrt(ns))
semeF_e = std(Ga_efah_e.individual/sqrt(ns))
semnR_e = std(Ga_nrah_e.individual/sqrt(ns))
semnF_e = std(Ga_nfah_e.individual/sqrt(ns))


timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000

r = figure;set(r,'position', [0 0 1000 400])
subplot(2,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erah_e.individual,1)), semeR_e, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efah_e.individual,1)), semeF_e, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time from amygdala early peaks (s)')
hold on
subplot(2,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrah_e.individual,1)), semnR_e, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfah_e.individual,1)), semnF_e, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time from amygdala early peaks (s)')
hold on
subplot(2,3,3)
plot(timestat,int_e.stat,'k','LineWidth',3)
xlabel('Time from amygdala early peaks (s)')
ylabel ('t-values')
%title('emotion by memory: EES early peaks amy hc')
hold on
subplot(2,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erah_l.individual,1)), semeR_l, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efah_l.individual,1)), semeF_l, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time from amygdala late peaks (s)')
hold on
subplot(2,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrah_l.individual,1)), semnR_l, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfah_l.individual,1)), semnF_l, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time from amygdala late peaks (s)')
hold on
subplot(2,3,6)
plot(timestat,int_l.stat,'k','LineWidth',3)
line([timestat(stat_line_l(1,1)), timestat(stat_line_l(1,2))],[-2, -2],'LineWidth',5, 'color', [0.9 0.5 0.6],'LineStyle','-'); 
xlabel('Time from amygdala late peaks (s)')
ylabel ('t-values')
