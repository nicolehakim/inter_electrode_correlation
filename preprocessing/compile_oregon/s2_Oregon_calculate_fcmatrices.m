clear a%% Analyze EEG connectivity from Oregon
clear all 
dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/';

load([dir,'data/preprocessed/oregon/Oregon_raw_compiled.mat'])
%% Settings 
nElec=size(data.erp.R_S2,3); 
nElec_anal=nElec-2;
nSubs=size(data.erp.R_S2,1);
nTrials=max([size(data.erp.R_S2,2),size(data.erp.L_C6,2),size(data.erp.R_C6,2),...
    size(data.erp.L_S2,2),size(data.erp.R_S6,2),size(data.erp.R_C2,2),...
    size(data.erp.L_S6,2),size(data.erp.L_C2,2)]);
Fs=data.srate; 
nTimes=size(data.erp.R_S2,4);
fc_mat.preTime=-200;%taken from Unsworth 2015
fc_mat.postTime=999;%taken from Unsworth 2015
fc_mat.sRate=250;%taken from Unsworth 2015
times=fc_mat.preTime:1000/fc_mat.sRate:fc_mat.postTime;%-data.preTime:data.postTime; 
baselineIntIdx=ismember(times,fc_mat.preTime:0);
retIntIdx=ismember(times,0:fc_mat.postTime);

nConds=4; 
nCondLabels={'C2','C6','S2','S6'};

tempEEGelec={'PO3', 'PO4', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'PO7', 'PO8', 'T7', 'T8', 'P7', 'P8', 'POz', 'Cz', 'Fz', 'Pz', 'HEOG', 'VEOG'};
data.tempEEGelec=tempEEGelec;
fc_mat.elecs_new=tempEEGelec(1:20);
fc_mat.elecs=data.chanLabels(1:20);
fc_mat.times=times;
%% First, make the organization of electrodes compatible with contra/ipsi design
% This basically involves flipping the order of the electrodes in one condition
% i.e. flip the electrode order for trials with stimuli presented on the right
% the same configuration as trials with stimuli presented on the left 

%pre-allocate some matrices
data.contraIpsiCorr.C2=nan(nSubs,nTrials*2,nElec,nTimes);%there will be some blank trials, but we can delete those at the end
data.contraIpsiCorr.C6=nan(nSubs,nTrials*2,nElec,nTimes);
data.contraIpsiCorr.S2=nan(nSubs,nTrials*2,nElec,nTimes);
data.contraIpsiCorr.S6=nan(nSubs,nTrials*2,nElec,nTimes);

%flip the electrodes in the left conditions
%Color - set size 2 condition 
data.contraIpsiCorr.C2(:,1:size(data.erp.L_C2,2),logical(mod(1:nElec,2)),:)=data.erp.L_C2(:,:,~mod(1:nElec,2),:,:);%flip the left to right
data.contraIpsiCorr.C2(:,1:size(data.erp.L_C2,2),logical(~mod(1:nElec,2)),:)=data.erp.L_C2(:,:,logical(mod(1:nElec,2)),:,:);%flip the right to left
data.contraIpsiCorr.C2(:,1:size(data.erp.L_C2,2),17:end,:)=data.erp.L_C2(:,:,17:end,:,:);%keep the non-lateralized ones the same
data.contraIpsiCorr.C2(:,size(data.erp.L_C2,2)+1:size(data.erp.L_C2,2)+size(data.erp.R_C2,2),:,:)=data.erp.R_C2;
data.contraIpsiCorr.C2=data.contraIpsiCorr.C2(:,~(data.contraIpsiCorr.C2(1,:,1,1)==0 | isnan(data.contraIpsiCorr.C2(1,:,1,1))),~ismember(tempEEGelec,{'HEOG','VEOG'}),:);%get rid of the excess trials & delete the extra electrodes (i.e. the horizontal and vertical eog channels)
%Color - set size 6 condition 
data.contraIpsiCorr.C6(:,1:size(data.erp.L_C6,2),logical(mod(1:nElec,2)),:)=data.erp.L_C6(:,:,~mod(1:nElec,2),:,:);%flip the left to right
data.contraIpsiCorr.C6(:,1:size(data.erp.L_C6,2),logical(~mod(1:nElec,2)),:)=data.erp.L_C6(:,:,logical(mod(1:nElec,2)),:,:);%flip the right to left
data.contraIpsiCorr.C6(:,1:size(data.erp.L_C6,2),17:end,:)=data.erp.L_C6(:,:,17:end,:,:);%keep the non-lateralized ones the same
data.contraIpsiCorr.C6(:,size(data.erp.L_C6,2)+1:size(data.erp.L_C6,2)+size(data.erp.R_C6,2),:,:)=data.erp.R_C6;
data.contraIpsiCorr.C6=data.contraIpsiCorr.C6(:,~(data.contraIpsiCorr.C6(1,:,1,1)==0 | isnan(data.contraIpsiCorr.C6(1,:,1,1))),~ismember(tempEEGelec,{'HEOG','VEOG'}),:);%get rid of the excess trials & delete the extra electrodes (i.e. the horizontal and vertical eog channels)
%Shape - set size 2 condition 
data.contraIpsiCorr.S2(:,1:size(data.erp.L_S2,2),logical(mod(1:nElec,2)),:)=data.erp.L_S2(:,:,~mod(1:nElec,2),:,:);%flip the left to right
data.contraIpsiCorr.S2(:,1:size(data.erp.L_S2,2),logical(~mod(1:nElec,2)),:)=data.erp.L_S2(:,:,logical(mod(1:nElec,2)),:,:);%flip the right to left
data.contraIpsiCorr.S2(:,1:size(data.erp.L_S2,2),17:end,:)=data.erp.L_S2(:,:,17:end,:,:);%keep the non-lateralized ones the same
data.contraIpsiCorr.S2(:,size(data.erp.L_S2,2)+1:size(data.erp.L_S2,2)+size(data.erp.R_S2,2),:,:)=data.erp.R_S2;
data.contraIpsiCorr.S2=data.contraIpsiCorr.S2(:,~(data.contraIpsiCorr.S2(1,:,1,1)==0 | isnan(data.contraIpsiCorr.S2(1,:,1,1))),~ismember(tempEEGelec,{'HEOG','VEOG'}),:);%get rid of the excess trials & delete the extra electrodes (i.e. the horizontal and vertical eog channels)
%Shape - set size 6 condition 
data.contraIpsiCorr.S6(:,1:size(data.erp.L_S6,2),logical(mod(1:nElec,2)),:)=data.erp.L_S6(:,:,~mod(1:nElec,2),:,:);%flip the left to right
data.contraIpsiCorr.S6(:,1:size(data.erp.L_S6,2),logical(~mod(1:nElec,2)),:)=data.erp.L_S6(:,:,logical(mod(1:nElec,2)),:,:);%flip the right to left
data.contraIpsiCorr.S6(:,1:size(data.erp.L_S6,2),17:end,:)=data.erp.L_S6(:,:,17:end,:,:);%keep the non-lateralized ones the same
data.contraIpsiCorr.S6(:,size(data.erp.L_S6,2)+1:size(data.erp.L_S6,2)+size(data.erp.R_S6,2),:,:)=data.erp.R_S6;
data.contraIpsiCorr.S6=data.contraIpsiCorr.S6(:,~(data.contraIpsiCorr.S6(1,:,1,1)==0 | isnan(data.contraIpsiCorr.S6(1,:,1,1))),~ismember(tempEEGelec,{'HEOG','VEOG'}),:);%get rid of the excess trials & delete the extra electrodes (i.e. the horizontal and vertical eog channels)

data.contraIpsiCorr.allTrial=cat(2,data.contraIpsiCorr.C2,data.contraIpsiCorr.C6,...
       data.contraIpsiCorr.S2,data.contraIpsiCorr.S6);%concatenate all of the data into one matrix
data.contraIpsiCorr.idx=[repmat(1,size(data.contraIpsiCorr.C2,2),1);...
    repmat(2,size(data.contraIpsiCorr.C6,2),1);...
    repmat(3,size(data.contraIpsiCorr.S2,2),1);...
    repmat(4,size(data.contraIpsiCorr.S6,2),1)];%make a condition index 

%% Low pass filter the eeg data at 50 hz to get rid of the noise from the crt monitor
tmp_contraIpsiCorr_C2=nan(size(data.contraIpsiCorr.C2));
tmp_contraIpsiCorr_C6=nan(size(data.contraIpsiCorr.C6));
tmp_contraIpsiCorr_S2=nan(size(data.contraIpsiCorr.S2));
tmp_contraIpsiCorr_S6=nan(size(data.contraIpsiCorr.S6));
tmp_contraIpsiCorr_allTrial=nan(size(data.contraIpsiCorr.allTrial));

tmp_erp_L_C2=nan(size(data.erp.L_C2));
tmp_erp_L_S6=nan(size(data.erp.L_S6));
tmp_erp_R_C2=nan(size(data.erp.R_C2));
tmp_erp_R_S6=nan(size(data.erp.R_S6));
tmp_erp_L_S2=nan(size(data.erp.L_S2));
tmp_erp_R_C6=nan(size(data.erp.R_C6));
tmp_erp_L_C6=nan(size(data.erp.L_C6));
tmp_erp_R_S2=nan(size(data.erp.R_S2));


for s=1:nSubs
    fprintf(['Filtering subject ',num2str(s),' out of ',num2str(nSubs),'\n'])
    for e=1:size(data.contraIpsiCorr.C2,3)
        for t=1:size(data.contraIpsiCorr.C2,2)
            tmp_contraIpsiCorr_C2(s,t,e,:)=eegfilt(squeeze(data.contraIpsiCorr.C2(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.contraIpsiCorr.C6,2)
            tmp_contraIpsiCorr_C6(s,t,e,:)=eegfilt(squeeze(data.contraIpsiCorr.C6(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.contraIpsiCorr.S2,2)
            tmp_contraIpsiCorr_S2(s,t,e,:)=eegfilt(squeeze(data.contraIpsiCorr.S2(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.contraIpsiCorr.S6,2)
            tmp_contraIpsiCorr_S6(s,t,e,:)=eegfilt(squeeze(data.contraIpsiCorr.S6(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.contraIpsiCorr.allTrial,2)
            tmp_contraIpsiCorr_allTrial(s,t,e,:)=eegfilt(squeeze(data.contraIpsiCorr.allTrial(s,t,e,:))',500,[],30);
        end
        
        for t=1:size(data.erp.L_C2,2)
            tmp_erp_L_C2(s,t,e,:)=eegfilt(squeeze(data.erp.L_C2(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.erp.L_S6,2)
            tmp_erp_L_S6(s,t,e,:)=eegfilt(squeeze(data.erp.L_S6(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.erp.R_C2,2)
            tmp_erp_R_C2(s,t,e,:)=eegfilt(squeeze(data.erp.R_C2(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.erp.R_S6,2)
            tmp_erp_R_S6(s,t,e,:)=eegfilt(squeeze(data.erp.R_S6(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.erp.L_S2,2)
            tmp_erp_L_S2(s,t,e,:)=eegfilt(squeeze(data.erp.L_S2(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.erp.R_C6,2)
            tmp_erp_R_C6(s,t,e,:)=eegfilt(squeeze(data.erp.R_C6(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.erp.L_C6,2)
            tmp_erp_L_C6(s,t,e,:)=eegfilt(squeeze(data.erp.L_C6(s,t,e,:))',500,[],30);
        end
        for t=1:size(data.erp.R_S6,2)
            tmp_erp_R_S2(s,t,e,:)=eegfilt(squeeze(data.erp.R_S6(s,t,e,:))',500,[],30);
        end
    end
end

% now, overwrite the original data, so that we only have the filtered data
data.contraIpsiCorr.C2=tmp_contraIpsiCorr_C2;
data.contraIpsiCorr.C6=tmp_contraIpsiCorr_C6;
data.contraIpsiCorr.S2=tmp_contraIpsiCorr_S2;
data.contraIpsiCorr.S6=tmp_contraIpsiCorr_S6;
data.contraIpsiCorr.allTrial=tmp_contraIpsiCorr_allTrial;

data.erp.L_C2=tmp_erp_L_C2;
data.erp.L_S6=tmp_erp_L_S6;
data.erp.R_C2=tmp_erp_R_C2;
data.erp.R_S6=tmp_erp_R_S6;
data.erp.L_S2=tmp_erp_L_S2;
data.erp.R_C6=tmp_erp_R_C6;
data.erp.L_C6=tmp_erp_L_C6;
data.erp.R_S2=tmp_erp_R_S2;
 
%% Now, do the pairwise pearson correlation separately for the color and shape conditions
% I do this to see whether we can predict performance on one task based on
% the eeg data from the other task 
fc_mat.baseline.shape=nan(nSubs,nElec_anal,nElec_anal); 
fc_mat.baseline.color=nan(nSubs,nElec_anal,nElec_anal); 
fc_mat.baseline.bothTasks=nan(nSubs,nElec_anal,nElec_anal); 

for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    
    tempColor=cat(1,squeeze(data.contraIpsiCorr.C2(s,:,:,baselineIntIdx)),squeeze(data.contraIpsiCorr.C6(s,:,:,baselineIntIdx)));
    tempColor=squeeze(nanmean(tempColor,1));
    tempShape=cat(1,squeeze(data.contraIpsiCorr.S2(s,:,:,baselineIntIdx)),squeeze(data.contraIpsiCorr.S6(s,:,:,baselineIntIdx)));
    tempShape=squeeze(nanmean(tempShape,1)); 
    tempD=squeeze(nanmean(data.contraIpsiCorr.allTrial(s,:,:,baselineIntIdx)));%average data over all trials to get: electrodes x time points
    %get connectivity matrix across subjects
    for e=1:nElec_anal
        for ee=1:nElec_anal
            if e==ee
                fc_mat.baseline.shape(s,e,ee)=19;
                fc_mat.baseline.color(s,e,ee)=19;
                fc_mat.baseline.bothTasks(s,e,ee)=19;
            else
                tempCol1=tempColor(e,:);
                tempCol2=tempColor(ee,:);
                tempShape1=tempShape(e,:); 
                tempShape2=tempShape(ee,:); 
                
                fc_mat.baseline.color(s,e,ee)=atanh(corr(tempCol1(:),tempCol2(:),'type','Pearson'));%correlate the two electrodes over time & Fisher normalize
                fc_mat.baseline.shape(s,e,ee)=atanh(corr(tempShape1(:),tempShape2(:),'type','Pearson'));
                
                temp1=tempD(e,:);
                temp2=tempD(ee,:);
                fc_mat.baseline.bothTasks(s,e,ee)=atanh(corr(temp1(:),temp2(:),'type','Pearson'));%correlate the two electrodes over time & Fisher normalize
            end
        end
    end
end

%% Now, do the pairwise pearson correlation separately for the color and shape conditions
% I do this to see whether we can predict performance on one task based on
% the eeg data from the other task 
fc_mat.RI.shape=nan(nSubs,nElec_anal,nElec_anal); 
fc_mat.RI.color=nan(nSubs,nElec_anal,nElec_anal); 
fc_mat.RI.bothTasks=nan(nSubs,nElec_anal,nElec_anal); 

for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    
    tempColor=cat(1,squeeze(data.contraIpsiCorr.C2(s,:,:,retIntIdx)),squeeze(data.contraIpsiCorr.C6(s,:,:,retIntIdx)));
    tempColor=squeeze(nanmean(tempColor,1));
    tempShape=cat(1,squeeze(data.contraIpsiCorr.S2(s,:,:,retIntIdx)),squeeze(data.contraIpsiCorr.S6(s,:,:,retIntIdx)));
    tempShape=squeeze(nanmean(tempShape,1)); 
    tempD=squeeze(nanmean(data.contraIpsiCorr.allTrial(s,:,:,retIntIdx)));%average data over all trials to get: electrodes x time points
    %get connectivity matrix across subjects
    for e=1:nElec_anal
        for ee=1:nElec_anal
            if e==ee
                fc_mat.RI.shape(s,e,ee)=19;
                fc_mat.RI.color(s,e,ee)=19;
                fc_mat.RI.bothTasks(s,e,ee)=19;
            else
                tempCol1=tempColor(e,:);
                tempCol2=tempColor(ee,:);
                tempShape1=tempShape(e,:); 
                tempShape2=tempShape(ee,:); 
                
                fc_mat.RI.color(s,e,ee)=atanh(corr(tempCol1(:),tempCol2(:),'type','Pearson'));%correlate the two electrodes over time & Fisher normalize
                fc_mat.RI.shape(s,e,ee)=atanh(corr(tempShape1(:),tempShape2(:),'type','Pearson'));
                
                temp1=tempD(e,:);
                temp2=tempD(ee,:);
                fc_mat.RI.bothTasks(s,e,ee)=atanh(corr(temp1(:),temp2(:),'type','Pearson'));%correlate the two electrodes over time & Fisher normalize
            end
        end
    end
end
%% Load the B2_int_closeavioral data 
% Not all of the eeg subjects have corresponding behavioral data, so we
% have to subset the data & only include the subjects with both eeg and
% behavioral data 

%load the behavioral data 
beh=load('/Volumes/eds-lab/_LongtermDataBackup/Kirsten/Data_BigStudy/BIG STUDY/OSF_ChangeDetection/K_cda_changedetection.mat');%these are the behavioral results from the eeg task
behidx=ismember(beh.cd_cda.subNum,data.subs);

fc_mat.k.ave=beh.cd_cda.K_ave(behidx);                      %average of all conditions k
fc_mat.k.color=nanmean(beh.cd_cda.K_cond(behidx,1:2),2);    %color k
fc_mat.k.shape=nanmean(beh.cd_cda.K_cond(behidx,3:4),2);    %shape k
fc_mat.k.color_ss2=beh.cd_cda.K_cond(behidx,1);             %color k, ss2
fc_mat.k.color_ss6=beh.cd_cda.K_cond(behidx,2);             %color k, ss6
fc_mat.k.shape_ss2=beh.cd_cda.K_cond(behidx,3);             %shape k, ss2
fc_mat.k.shape_ss6=beh.cd_cda.K_cond(behidx,4);             %shape k, ss6

%% downsample the fc_mat, so that we only have subjects that have both eeg and behavioral data
eegidx=ismember(data.subs,beh.cd_cda.subNum);
fc_mat.baseline.shape=fc_mat.baseline.shape(eegidx,:,:);
fc_mat.baseline.color=fc_mat.baseline.color(eegidx,:,:);
fc_mat.baseline.bothTasks=fc_mat.baseline.bothTasks(eegidx,:,:); 

fc_mat.RI.shape=fc_mat.RI.shape(eegidx,:,:);
fc_mat.RI.color=fc_mat.RI.color(eegidx,:,:);
fc_mat.RI.bothTasks=fc_mat.RI.bothTasks(eegidx,:,:); 

%% Save the correlation matrix
save([dir,'/data/compiled/','oregon_site_fc_mat.mat'],'fc_mat','-v7.3');
% save([path,'BigStudyData_compiled.mat'],'data','-v7.3');%save the raw data structure again because we did the contra/ipsi transformation 

