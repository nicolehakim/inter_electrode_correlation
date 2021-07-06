%% Compile the data from the big study to use in the CPM analysis
%
% Written by Nicole Hakim, November 2019 
%% 
clear all 
%% Get some basic information 
dataPath='/Volumes/eds-lab/_LongtermDataBackup/Kirsten/Data_BigStudy/BIG STUDY/CDA_data/';%location of the raw data

dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/';%
savePath=[dir,'/data/preprocessed/oregon/'];
subs=[1:8,10:18,20:26,28:34,36:47,49,51:52,54:59,61:63,65:68,70,73,74,78:102,...
    104:112,115:121,123:133,135:138,140:147,150,152:156,158:160,162:171,...
    174:175,178:180,184:191,194,197,198,203,206,208:219];
nSubs=length(subs);
data.subs=subs;
nElec=22; 
nTimes=300;
nTrialsMax=100;
fileName='/erp_singletrial.mat';

%% Get the data for each individual participant
data.erp.L_C2=nan(nSubs,nTrialsMax,nElec,nTimes);
data.erp.L_S6=nan(nSubs,nTrialsMax,nElec,nTimes);
data.erp.R_C2=nan(nSubs,nTrialsMax,nElec,nTimes);
data.erp.R_S6=nan(nSubs,nTrialsMax,nElec,nTimes);
data.erp.L_S2=nan(nSubs,nTrialsMax,nElec,nTimes);
data.erp.R_C6=nan(nSubs,nTrialsMax,nElec,nTimes);
data.erp.L_C6=nan(nSubs,nTrialsMax,nElec,nTimes);
data.erp.R_S2=nan(nSubs,nTrialsMax,nElec,nTimes);

for s=1:nSubs
    fprintf(['Subject ',num2str(s),' out of ', num2str(nSubs),'\n'])
    sn=subs(s);
    load([dataPath,num2str(sn),fileName]);
    
    data.erp.L_C2(s,1:size(erp.trial.L_C2,1),:,:)=erp.trial.L_C2; 
    data.erp.L_S6(s,1:size(erp.trial.L_S6,1),:,:)=erp.trial.L_S6; 
    data.erp.R_C2(s,1:size(erp.trial.R_C2,1),:,:)=erp.trial.R_C2; 
    data.erp.R_S6(s,1:size(erp.trial.R_S6,1),:,:)=erp.trial.R_S6; 
    data.erp.L_S2(s,1:size(erp.trial.L_S2,1),:,:)=erp.trial.L_S2; 
    data.erp.R_C6(s,1:size(erp.trial.R_C6,1),:,:)=erp.trial.R_C6; 
    data.erp.L_C6(s,1:size(erp.trial.L_C6,1),:,:)=erp.trial.L_C6; 
    data.erp.R_S2(s,1:size(erp.trial.R_S2,1),:,:)=erp.trial.R_S2; 
    
end

data.chanLabels=erp.allChans; 
data.preTime=erp.pre_timepoint;
data.postTime=erp.post_timepoint;
data.srate=erp.srate;

%% Save
save([savePath,'Oregon_raw_compiled.mat'],'data','-v7.3')

