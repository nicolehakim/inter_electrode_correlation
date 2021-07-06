%% This trains and tests a model on the non-flipped Oregon data 
clear all
%% Settings
dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/';%where the data is saved
nFolds=5; 
nIter=10000;

%% Load the data
load([dir,'data/compiled/oregon_site_fc_mat_nonflipped_electrodes.mat']);

fc_both=bs.fc_mat.overlappingElecs.RI.bothTasks;
fc_beh.both=bs.fc_mat.k.ave;

%Need to realign electrodes and only include relevant ones
nSubs=size(fc_both,1);
subs=1:nSubs; 
nElec=size(fc_both,2); 

%Which edges do you want to consider to be important? 
perc_importEdges=.10;% This is the percentage of top edges that will be included (i.e. if perc_importEdges=20, then 20% of the top edges will be included 
total_importantEdges=round(((nElec*nElec-nElec)*perc_importEdges)/2);%number of positive and number of negative edges to include

nSubs_perFold=floor(nSubs/nFolds);
%% Train and test the model
z_across=nan(nIter,1);
mse_across=nan(nIter,1);

z_shuff=nan(nIter,1);
mse_shuff=nan(nIter,1);

y_pre=nan(nIter,nFolds,nSubs_perFold);
y_act=nan(nIter,nFolds,nSubs_perFold);

for ii =1:nIter
    fprintf(['Iteration ',num2str(ii),' out of ',num2str(nIter),'\n'])
    
    subIdx=shuffle(subs); 
    subs_temp=subIdx(1:nSubs_perFold*nFolds);
    foldIdx=reshape(subs_temp(1:(nSubs_perFold*nFolds)),nFolds,nSubs_perFold);
    
    z_temp=nan(nFolds,1); 
    mse_temp=nan(nFolds,1);
    
    z_temp_shuff=nan(nFolds,1); 
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
        z_temp(k)=atanh(corr(y_predict,y_actual,'type','Spearman'));
        mse_temp(k)=immse(y_predict,y_actual);
        %shuffle
        shuff_actual=shuffle(y_actual); 
        z_temp_shuff(k)=atanh(corr(y_predict,shuff_actual,'type','Spearman'));
        mse_temp_shuff(k)=immse(y_predict,shuff_actual);
    end
    
    %corr
    z_across(ii)=nanmean(z_temp); 
    mse_across(ii)=nanmean(mse_temp); 
    %shuffle
    z_shuff(ii)=nanmean(z_temp_shuff); 
    mse_shuff(ii)=nanmean(mse_temp_shuff); 
end

%% calculate non-paramemtric p-values
r_across=tanh(z_across);
mean_r=tanh(mean(z_across));
median_r=tanh(median(z_across));
r_shuff=tanh(z_shuff);

p_para=(1+sum(r_shuff>=mean_r))/(1+nIter);
mse_ave=nanmean(mse_across);
p_para_mse=(1+sum(mse_shuff<=mean(mse_across)))/(1+nIter);

%% Finished, and let me know
load handel
sound(y,Fs)
