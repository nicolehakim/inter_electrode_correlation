%% Big study analysis: K fold cross validation
% this script does K-fold CPM only with the averaged fc matrix
% also runs the K-fold CPM with raw eeg (not connectivity)
%% Train and test CPM with connectivity matrix data 
clear all 
%% Paths 
dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/';
%% Load the data 
%% Settings to change
nFolds=5; 
nIter=10000; 

load([dir,'/data/compiled/oregon_site_fc_mat.mat'])%load the data
fc_color=bs.fc_mat.overlappingElecs.RI.color;
fc_shape=bs.fc_mat.overlappingElecs.RI.shape;

nElec=size(fc_color,2);

perc_importEdges=10;% This is the percentage of top edges that will be included (i.e. if perc_importEdges=10, then 5% of the top positive and 5% of the top negative edges will be included 
total_importantEdges=round(((nElec*nElec-nElec)*perc_importEdges)/2);%number of positive and number of negative edges to include

%% Settings that don't change
fc_beh.color=bs.fc_mat.k.color;
fc_beh.shape=bs.fc_mat.k.shape;
fc_beh.both=bs.fc_mat.k.ave;

nSubs=size(fc_color,1);
subs=1:nSubs;

nSubs_perFold=floor(nSubs/nFolds);
%% Train and test CPM (K fold) - on functional connectivity data (averaged over the RI times) - with nonparametric test of value 
y_predict.trainCol=nan(nIter,nFolds,nSubs_perFold); 
y_actual.trainCol=nan(nIter,nFolds,nSubs_perFold,1);
y_predict.trainColTestShape=nan(nIter,nFolds,nSubs_perFold); 

y_predict.trainShape=nan(nIter,nFolds,nSubs_perFold); %first dimension: train shape, test shape. second dimension: train shape, test color
y_actual.trainShape=nan(nIter,nFolds,nSubs_perFold,1);
y_predict.trainShapeTestCol=nan(nIter,nFolds,nSubs_perFold);

y_subs=nan(nIter,nFolds,nSubs_perFold);

for ii=1:nIter
    fprintf(['Iteration ',num2str(ii),' out of ',num2str(nIter),'\n'])
    subIdx=shuffle(subs); foldIdx=reshape(subIdx(1:nFolds*nSubs_perFold),[nFolds,nSubs_perFold]);
    for f=1:nFolds
        fprintf(['Fold ',num2str(f),' out of ',num2str(nFolds),'\n'])
        testSub=foldIdx(f,:);
        trainSubs=foldIdx(foldIdx~=testSub);
        
        rhoEdges=nan(2,nElec,nElec); 
        for ed1=1:nElec
            for ed2=1:nElec
                %train color
                [r,p]=corr(squeeze(fc_color(trainSubs,ed1,ed2)),fc_beh.color(trainSubs));
                rhoEdges(1,ed1,ed2)=r;
                
                %train shape
                [r,p]=corr(squeeze(fc_shape(trainSubs,ed1,ed2)),fc_beh.shape(trainSubs));
                rhoEdges(2,ed1,ed2)=r;
            end
        end

        importantEdges=nan(2,nElec,nElec);
        
        %color task
        tmp=rhoEdges(1,:,:);
        tmp=sort(tmp(~isnan(tmp)),1,'ascend');
        importantEdges(1,squeeze(ismember(rhoEdges(1,:,:),tmp(1:perc_importEdges))))=-1;%negative edges
        importantEdges(1,squeeze(ismember(rhoEdges(1,:,:),tmp(end-perc_importEdges:end))))=1;%positive edges
        
        %shape task
        tmp=rhoEdges(2,:,:);
        tmp=sort(tmp(~isnan(tmp)),1,'ascend');
        importantEdges(2,squeeze(ismember(rhoEdges(2,:,:),tmp(1:perc_importEdges))))=-1;%negative edges
        importantEdges(2,squeeze(ismember(rhoEdges(2,:,:),tmp(end-perc_importEdges:end))))=1;%positive edges
        
        %for each subject, summarize selected features
        summary_features=nan(nSubs,2); %1st dimension is train color, test shape; 2nd dimension is train shape, test color
        for i=1:nSubs
            thisSub_color=squeeze(fc_color(i,:,:));
            summary_features(i,1)=nanmean(thisSub_color(squeeze(importantEdges(1,:,:))>0))-nanmean(thisSub_color(squeeze(importantEdges(1,:,:))<0));
            
            thisSub_shape=squeeze(fc_shape(i,:,:));
            summary_features(i,2)=nanmean(thisSub_shape(squeeze(importantEdges(2,:,:))>0))-nanmean(thisSub_shape(squeeze(importantEdges(2,:,:))<0));
        end
        
        if isnan(mean(summary_features(:,1),1)) | isnan(mean(summary_features(:,2),1))%if there are not significant positive or negative edges in this iteration, then skip the rest of this iteration's loop 
            continue
        end
        
        %train across subjects for each task separately
        mdl_col=robustfit(summary_features(trainSubs,1),fc_beh.color(trainSubs));
        mdl_shape=robustfit(summary_features(trainSubs,2),fc_beh.shape(trainSubs));
        
        %test held out subjects
        y_predict.trainCol(ii,f,:)=mdl_col(2)*summary_features(testSub,1)+mdl_col(1);%train color, test color
        y_actual.trainCol(ii,f,:)=fc_beh.color(testSub);
        y_predict.trainColTestShape(ii,f,:)=mdl_col(2)*summary_features(testSub,2)+mdl_col(1);%train color, test shape
        
        y_predict.trainShape(ii,f,:)=mdl_shape(2)*summary_features(testSub,2)+mdl_shape(1);
        y_actual.trainShape(ii,f,:)=fc_beh.shape(testSub);
        y_predict.trainShapeTestCol(ii,f,:)=mdl_shape(2)*summary_features(testSub,1)+mdl_shape(1);%train shape, test color
        
        y_subs(ii,f,:)=testSub;
        
        
    end
end
%% Average the prediction over iterations for each subject
z_across_trainC_testS=nan(nIter,1);
z_shuff_trainC_testS=nan(nIter,1);
z_across_trainS_testC=nan(nIter,1);
z_shuff_trainS_testC=nan(nIter,1);

mse_across_trainC_testS=nan(nIter,1);
mse_shuff_trainC_testS=nan(nIter,1);
mse_across_trainS_testC=nan(nIter,1);
mse_shuff_trainS_testC=nan(nIter,1);

for ii=1:nIter
    fprintf(['iteration ',num2str(ii),' out of ',num2str(nIter),'\n'])
    z_temp_trainC_testS=nan(1,nFolds);
    z_temp_shuff_trainC_testS=nan(1,nFolds);
    mse_temp_trainC_testS=nan(1,nFolds);
    mse_temp_shuff_trainC_testS=nan(1,nFolds);
    
    z_temp_trains_testc=nan(1,nFolds);
    z_temp_shuff_trains_testc=nan(1,nFolds);
    mse_temp_trains_testc=nan(1,nFolds);
    mse_temp_shuff_trains_testc=nan(1,nFolds);
    
    for f=1:nFolds
        z_temp_trainC_testS(f)=atanh(corr(squeeze(y_predict.trainColTestShape(ii,f,:)),squeeze(y_actual.trainShape(ii,f,:)),'type','Spearman'));
        z_temp_shuff_trainC_testS(f)=atanh(corr(shuffle(squeeze(y_predict.trainColTestShape(ii,f,:))),squeeze(y_actual.trainShape(ii,f,:)),'type','Spearman'));
        mse_temp_trainC_testS(f)=immse(squeeze(y_predict.trainColTestShape(ii,f,:)),squeeze(y_actual.trainShape(ii,f,:)));
        mse_temp_shuff_trainC_testS(f)=immse(shuffle(squeeze(y_predict.trainColTestShape(ii,f,:))),squeeze(y_actual.trainShape(ii,f,:)));
        
        z_temp_trains_testc(f)=atanh(corr(squeeze(y_predict.trainShapeTestCol(ii,f,:)),squeeze(y_actual.trainCol(ii,f,:)),'type','Spearman'));
        z_temp_shuff_trains_testc(f)=atanh(corr(shuffle(squeeze(y_predict.trainShapeTestCol(ii,f,:))),squeeze(y_actual.trainCol(ii,f,:)),'type','Spearman'));
        mse_temp_trains_testc(f)=immse(squeeze(y_predict.trainShapeTestCol(ii,f,:)),squeeze(y_actual.trainCol(ii,f,:)));
        mse_temp_shuff_trains_testc(f)=immse(shuffle(squeeze(y_predict.trainShapeTestCol(ii,f,:))),squeeze(y_actual.trainCol(ii,f,:)));
    end
    
    z_across_trainC_testS(ii)=nanmean(z_temp_trainC_testS);
    z_shuff_trainC_testS(ii)=nanmean(z_temp_shuff_trainC_testS);
    mse_across_trainC_testS(ii)=nanmean(mse_temp_trainC_testS);
    mse_shuff_trainC_testS(ii)=nanmean(mse_temp_shuff_trainC_testS);
    
    
    z_across_trainS_testC(ii)=nanmean(z_temp_trains_testc);
    z_shuff_trainS_testC(ii)=nanmean(z_temp_shuff_trains_testc);
    mse_across_trainS_testC(ii)=nanmean(mse_temp_trains_testc);
    mse_shuff_trainS_testC(ii)=nanmean(mse_temp_shuff_trains_testc);
end


%% calculate non-paramemtric p-values: train color, test shape
r_across_trainC_testS=tanh(z_across_trainC_testS);
mean_r_trainC_testS=tanh(mean(z_across_trainC_testS));
median_r_trainC_testS=tanh(median(z_across_trainC_testS));
r_shuff_trainC_testS=tanh(z_shuff_trainC_testS);

p_para_trainC_testS=(1+sum(r_shuff_trainC_testS>=mean_r_trainC_testS))/(1+nIter);
mse_ave_trainC_testS=nanmean(mse_across_trainC_testS);
p_para_mse_trainC_testS=(1+sum(mse_shuff_trainC_testS<=mean(mse_across_trainC_testS)))/(1+nIter);

%% calculate non-paramemtric p-values: train shape, test color
r_across_trainS_testC=tanh(z_across_trainS_testC);
mean_r_trainS_testC=tanh(mean(z_across_trainS_testC));
median_r_trainS_testC=tanh(median(z_across_trainS_testC));
r_shuff_trainS_testC=tanh(z_shuff_trainS_testC);

p_para_trainS_testC=(1+sum(r_shuff_trainS_testC>=mean_r_trainS_testC))/(1+nIter);
mse_ave_trainS_testC=nanmean(mse_across_trainS_testC);
p_para_mse_trainS_testC=(1+sum(mse_shuff_trainS_testC<=mean(mse_across_trainS_testC)))/(1+nIter);
