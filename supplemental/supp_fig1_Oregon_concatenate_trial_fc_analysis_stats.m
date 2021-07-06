%% This trains and tests a model on the concatenated data for the Oregon data set 
clear all 
%% Load the data
nFolds=5; 
nIter=10000;

dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg';

load([dir,'/data/compiled/oregon_concatenated_trials.mat']);

fc_both=bs.fc_mat.overlappingElecs.RI.bothTasks;
fc_beh.both=bs.fc_mat.k.ave;

%Need to realign electrodes and only include relevant ones
nSubs=size(fc_both,1);
subs=1:nSubs; 
nElec=size(fc_both,2); 

%Which edges do you want to consider to be important? 
perc_importEdges=.10;% This is the percentage of top edges that will be included (i.e. if perc_importEdges=20, then 20% of the top edges will be included 
total_importantEdges=round(((nElec*nElec-nElec)*perc_importEdges)/2);%number of positive and number of negative edges to include
%% Train and test the model
r_across=nan(nIter,1);
mse_across=nan(nIter,1);

r_shuff=nan(nIter,1);
mse_shuff=nan(nIter,1);

y_pre=nan(nIter,nFolds,floor(nSubs/nFolds));
y_act=nan(nIter,nFolds,floor(nSubs/nFolds));

for ii =1:nIter
    fprintf(['Iteration ',num2str(ii),' out of ',num2str(nIter),'\n'])
    
    subIdx=shuffle(subs); 
    subs_temp=subIdx(1:floor(nSubs/nFolds)*nFolds);
    foldIdx=reshape(subs_temp(1:(floor(nSubs/nFolds)*nFolds)),nFolds,floor(nSubs/nFolds));
    
    r_temp=nan(nFolds,1); 
    mse_temp=nan(nFolds,1);
    
    r_temp_shuff=nan(nFolds,1); 
    mse_temp_shuff=nan(nFolds,1);
    
    for k=1:nFolds
        %fprintf(['Fold ',num2str(k),' out of ',num2str(nFolds),'\n'])
        
        testSub=foldIdx(k,:);
        trainSubs=subs(~ismember(subs,testSub));
        
        rhoEdges=nan(nElec,nElec);
        for ed1=1:nElec
            for ed2=1:nElec
                [r,p]=corr(squeeze(fc_both(trainSubs,ed1,ed2)),fc_beh.both(trainSubs),'type','Spearman');
                rhoEdges(ed1,ed2)=r;
            end
        end
        
        %Get the important edges
        importantEdges=nan(nElec,nElec);
        %color task
        tmp=squeeze(rhoEdges(:,:));
        tmp=sort(tmp(~isnan(tmp)),1,'ascend');
        importantEdges(squeeze(ismember(rhoEdges(:,:),tmp(1:total_importantEdges))))=-1;%negative edges
        importantEdges(squeeze(ismember(rhoEdges(:,:),tmp(end-total_importantEdges:end))))=1;%positive edges
        
        %for each subject, summarize selected features
        summary_features=nan(nSubs,1); %1st dimension is train color, test shape; 2nd dimension is train shape, test color
        for i=1:nSubs
            thisSub=squeeze(fc_both(i,:,:));
            summary_features(i,1)=nanmean(thisSub(squeeze(importantEdges(:,:))>0))-nanmean(thisSub(squeeze(importantEdges(:,:))<0));
        end
       
        %train across subjects for each task separately
        mdl=robustfit(summary_features(trainSubs,1),fc_beh.both(trainSubs));
        
        %test held out subjects
        y_actual=fc_beh.both(testSub);
        %make prediction
        y_predict=mdl(2)*summary_features(testSub,1)+mdl(1);
        
        %save this for later
        y_pre(ii,k,:)=y_predict;
        y_act(ii,k,:)=y_actual;
        
        %calculate results 
        r_temp(k)=corr(y_predict,y_actual,'type','Pearson');
        mse_temp(k)=immse(y_predict,y_actual);
        %shuffle
        shuff_actual=shuffle(y_actual); 
        r_temp_shuff(k)=corr(y_predict,shuff_actual,'type','Pearson');
        mse_temp_shuff(k)=immse(y_predict,shuff_actual);
    end
    
    %corr
    r_across(ii)=nanmean(r_temp); 
    mse_across(ii)=nanmean(mse_temp); 
    %shuffle
    r_shuff(ii)=nanmean(r_temp_shuff); 
    mse_shuff(ii)=nanmean(mse_temp_shuff); 
end
%% calculate non-paramemtric p-value
p_para=(1+sum(r_shuff>=r_across))/(1+nIter);
cohens_d=computeCohen_d(r_shuff,r_across,'paired');
mean_r=nanmean(r_across); 
mean_mse=nanmean(mse_across);