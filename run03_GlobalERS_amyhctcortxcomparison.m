
clear all
close all
clc

restoredefaultpath
addpath '~/fieldtrip-20210212/'
ft_defaults

%% reorganize tmp cortex electrodes order as amy and hc ones
subjectsamyhc = [60 13 15 160 25 32 33 600 2000 36 38 8]; 
subjectscrtx=   [60 13 15 160 25 32 33 36 38 2000 600 8];

cd ~/data2github
load('Ga_GlobalERSCrtx.mat')

tmper=Ga_er
tmpekf=Ga_ekf
tmpnr=Ga_nr
tmpnkf=Ga_nkf
%%

Ga_er.powspctrm([8:12],:,:,:)=[]
Ga_er.powspctrm(8,:,:,:)=tmper.powspctrm(11,:,:,:);
Ga_er.powspctrm(9,:,:,:)=tmper.powspctrm(10,:,:,:);
Ga_er.powspctrm(10,:,:,:)=tmper.powspctrm(8,:,:,:);
Ga_er.powspctrm(11,:,:,:)=tmper.powspctrm(9,:,:,:);
Ga_er.powspctrm(12,:,:,:)=tmper.powspctrm(12,:,:,:);

Ga_ekf.powspctrm([8:12],:,:,:)=[]
Ga_ekf.powspctrm(8,:,:,:)=tmpekf.powspctrm(11,:,:,:);
Ga_ekf.powspctrm(9,:,:,:)=tmpekf.powspctrm(10,:,:,:);
Ga_ekf.powspctrm(10,:,:,:)=tmpekf.powspctrm(8,:,:,:);
Ga_ekf.powspctrm(11,:,:,:)=tmpekf.powspctrm(9,:,:,:);
Ga_ekf.powspctrm(12,:,:,:)=tmpekf.powspctrm(12,:,:,:);

Ga_nr.powspctrm([8:12],:,:,:)=[]
Ga_nr.powspctrm(8,:,:,:)=tmpnr.powspctrm(11,:,:,:);
Ga_nr.powspctrm(9,:,:,:)=tmpnr.powspctrm(10,:,:,:);
Ga_nr.powspctrm(10,:,:,:)=tmpnr.powspctrm(8,:,:,:);
Ga_nr.powspctrm(11,:,:,:)=tmpnr.powspctrm(9,:,:,:);
Ga_nr.powspctrm(12,:,:,:)=tmpnr.powspctrm(12,:,:,:);

Ga_nkf.powspctrm([8:12],:,:,:)=[]
Ga_nkf.powspctrm(8,:,:,:)=tmpnkf.powspctrm(11,:,:,:);
Ga_nkf.powspctrm(9,:,:,:)=tmpnkf.powspctrm(10,:,:,:);
Ga_nkf.powspctrm(10,:,:,:)=tmpnkf.powspctrm(8,:,:,:);
Ga_nkf.powspctrm(11,:,:,:)=tmpnkf.powspctrm(9,:,:,:);
Ga_nkf.powspctrm(12,:,:,:)=tmpnkf.powspctrm(12,:,:,:);
%% select mat values based on t and f of the amy and hc effect

subjects = [60 13 15 160 25 32 33 600 2000 36 38 8]; 
amymat=load('amymat.mat')
hcmat=load('hcmat.mat')

toi=[min(amymat.tamy(1), hcmat.t(1)) max(amymat.tamy(2), hcmat.t(2))]
foi=[min(amymat.famy(1), hcmat.f(1)) max(amymat.famy(2), hcmat.f(2))]

t= toi
f = foi

pt1 = nearest(Ga_er.time,t(1));
pt2 = nearest(Ga_er.time,t(2));
pf1 = nearest(Ga_er.freq,f(1));
pf2 = nearest(Ga_er.freq,f(2));

%clear mat 
crtxmat(:,1) = squeeze(mean(mean(Ga_er.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
crtxmat(:,2) = squeeze(mean(mean(Ga_ekf.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
crtxmat(:,3) = squeeze(mean(mean(Ga_nr.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
crtxmat(:,4) = squeeze(mean(mean(Ga_nkf.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
%% compute the interaction in each ROI

inthc=[hcmat.mat(:,1)-hcmat.mat(:,2)]-[hcmat.mat(:,3)-hcmat.mat(:,4)];
intamy=[amymat.mat(:,1)-amymat.mat(:,2)]-[amymat.mat(:,3)-amymat.mat(:,4)];
intcrtx=[crtxmat(:,1)-crtxmat(:,2)]-[crtxmat(:,3)-crtxmat(:,4)];

all=[intamy,inthc,intcrtx]
%%
[Hah,Pah,CIah,STATSah] = ttest(all(:,1),all(:,2));
[Hhc,Phc,CIhc,STATShc] = ttest(all(:,2),all(:,3));
[Hac,Pac,CIac,STATSac] = ttest(all(:,1),all(:,3));
%% plot Supplementary Fig. 7
addpath('~/utils/beeswarm-master')
addpath('~/utils/stdshade.m')

figure
ns=size(all,1)
x = [ones(ns,1) ones(ns,1)*2 ones(ns,1)*3];
y = [all(:,1) all(:,2) all(:,3)];
beeswarm(x(:),y(:),'sort_style','up','dot_size',2,'overlay_style','sd','colormap',[1 0.4 0.7; 0 0 1;0 1 0])
%ylim([-0.2 0.2]);
ylabel('mean ERS values','FontSize',10)
xticklabels({[],'Amygdala',[],'Hippocampus',[],'LatCortex'})
%%
meanval = mean(all);
seval =std(all)/sqrt(size(subjects,2))