clear all
close all
clc

%%
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
toi=[-0.5 1.5]
%% find amygdala peaks during encoding and prepare the data for the amygdala encoding-retrieval similarity
parfor s=1:12;
[TFpeakall_era{s}, TFpeakall_erh{s}, TFpeakall_ehith{s}, nb_er{s}] = findamypeaksEESERS(eR_amy{s},TF_era(s),Amylist(s), TF_era(s),TF_ehita(s),Amylist(s),50,freq,toi);
[TFpeakall_nra{s}, TFpeakall_nrh{s}, TFpeakall_nhith{s}, nb_nr{s}] = findamypeaksEESERS(nR_amy{s},TF_nra(s),Amylist(s), TF_nra(s),TF_nhita(s),Amylist(s),50,freq,toi);
[TFpeakall_efa{s}, TFpeakall_efh{s}, TFpeakall_emissh{s}, nb_ef{s}] = findamypeaksEESERS(eKF_amy{s},TF_efa(s),Amylist(s), TF_efa(s),TF_emissa(s),Amylist(s),50,freq,toi);
[TFpeakall_nfa{s}, TFpeakall_nfh{s}, TFpeakall_nmissh{s}, nb_nf{s}] = findamypeaksEESERS(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfa(s),TF_nmissa(s),Amylist(s),50,freq,toi);
end

parfor s=1:12;
[TFpeakall_erac1{s}, TFpeakall_erhc1{s}, TFpeakall_ehithc1{s}, nb_er{s}] = findamypeaksEESERS(eR_amy{s},TF_era(s),Amylist(s), TF_era(s),TF_ehita(s),Amylist(s),150,freq,toi);
[TFpeakall_nrac1{s}, TFpeakall_nrhc1{s}, TFpeakall_nhithc1{s}, nb_nr{s}] = findamypeaksEESERS(nR_amy{s},TF_nra(s),Amylist(s), TF_nra(s),TF_nhita(s),Amylist(s),150,freq,toi);
[TFpeakall_efac1{s}, TFpeakall_efhc1{s}, TFpeakall_emisshc1{s}, nb_ef{s}] = findamypeaksEESERS(eKF_amy{s},TF_efa(s),Amylist(s), TF_efa(s),TF_emissa(s),Amylist(s),150,freq,toi);
[TFpeakall_nfac1{s}, TFpeakall_nfhc1{s}, TFpeakall_nmisshc1{s}, nb_nf{s}] = findamypeaksEESERS(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfa(s),TF_nmissa(s),Amylist(s),150,freq,toi);
end

%%
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_ehith{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]
%%
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

%%
load Gatemp.mat

GAs = GTFs_Multitaper_baseline_Unpl
GAs = rmfield(GAs,'powspctrm');
GAs = rmfield(GAs,'freq');


GAs.time = []
dtp= size(alleraht_ers{1},2)-1
GAs.time = TF_erh(1).values.time(151:151+dtp);

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:151+dtp);

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

%%
URUFahers = Ga_erahers;
URUFahers.individual = Ga_erahers.individual - Ga_efahers.individual
NRNFahers = Ga_erahers;
NRNFahers.individual = Ga_nrahers.individual - Ga_nfahers.individual

URUFahc1ers = Ga_erahc1ers;
URUFahc1ers.individual = Ga_erahc1ers.individual - Ga_efahc1ers.individual
NRNFahc1ers = Ga_erahc1ers;
NRNFahc1ers.individual = Ga_nrahc1ers.individual - Ga_nfahc1ers.individual

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

ns=size(URUFahers.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

int_ers = ft_timelockstatistics(cfg,URUFahers, NRNFahers)
intc1_ers= ft_timelockstatistics(cfg,URUFahc1ers, NRNFahc1ers)

%% to plot results in Supplementary Fig.8

cd ~/data2github
load('statERS_Amygdala.mat'); %from 0 to 1.5


ns=size(Ga_erahers.individual,1)
%% plot Supplementary Fig.8

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
            stat_linec1(i,:) = [wherec1(1) wherec1(end)];
end
%% 
ns=size(Ga_erahers.individual,1)

semeRers = std(Ga_erahers.individual/sqrt(ns))
semeFers = std(Ga_efahers.individual/sqrt(ns))
semnRers = std(Ga_nrahers.individual/sqrt(ns))
semnFers = std(Ga_nfahers.individual/sqrt(ns))

semeRc1 = std(Ga_erahc1ers.individual/sqrt(ns))
semeFc1 = std(Ga_efahc1ers.individual/sqrt(ns))
semnRc1 = std(Ga_nrahc1ers.individual/sqrt(ns))
semnFc1 = std(Ga_nfahc1ers.individual/sqrt(ns))


timestat=[0:21:1500];
timestat=timestat(1,[1:71])/1000
time=[-500:21:1500];
time=time(1,[1:96])/1000

r=figure;set(r,'position', [0 0 1000 400])
subplot(2,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahers.individual,1)), semeRers, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahers.individual,1)), semeFers, 'm', 1); hold on
%legend('erhit','ekhit&miss')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(2,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahers.individual,1)), semnRers, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahers.individual,1)), semnFers, 'k', 1); hold on
%legend('nrhit','nkhit&miss')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(2,3,3)
plot(timestat,int_ers.stat,'k','LineWidth',3)
hold on
xlabel('Time (s)')
ylabel ('t-values')
hold on
subplot(2,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erahc1ers.individual,1)), semeRc1, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahc1ers.individual,1)), semeFc1, 'm', 1); hold on
%legend('erhit','ekhit&miss')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(2,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrahc1ers.individual,1)), semnRc1, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahc1ers.individual,1)), semnFc1, 'k', 1); hold on
%legend('nrhit','nkhit&miss')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(2,3,6)
plot(timestat,intc1_ers.stat,'k','LineWidth',3)
hold on
xlabel('Time (s)')
ylabel ('t-values')

