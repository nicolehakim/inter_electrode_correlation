%% Train and test gF with 5-fold CV, using its own edges
clear all; close all; 
%% Load the data 
dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/';

% do these analyses with the edges that overlapped 
% load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/oregon_chicago_top10percedges.mat')%sig_edges_trainOreg_testChi.mat')
% importantEdges=sigEdges_overlap.oregon_chicago;%do the  only with the overlapping edges from Oregon and Chicago for both datasets
% 
% eeg_subs=data.subs(ismember(data.subs,beh.cd_cda.subNum));
% eeg_sub_idx=ismember(eeg_subs,dat.subjects); 
% behav_other_subs=ismember(dat.subjects,eeg_subs);

load([dir,'/data/compiled/chicago_site_fc_mat.mat'])

fc_both    = fc_mat.eeg;
behav = fc_mat.beh;

nElec  = size(fc_both,2);%17;   % number of nodes

nIter=10000; 
nFolds=5; 
subs=1:size(fc_both,1);
nSubs=length(subs); 

%Which edges do you want to consider to be important? 
perc_importEdges=.10;% This is the percentage of top edges that will be included (i.e. if perc_importEdges=20, then 20% of the top edges will be included 
total_importantEdges=round(((nElec*nElec-nElec)*perc_importEdges)/2);%number of positive and number of negative edges to include

%% Predict gF
z_across=nan(nIter,1);
mse_across=nan(nIter,1);

z_shuff=nan(nIter,1);
mse_shuff=nan(nIter,1);

y_pre=nan(nIter,nFolds,floor(nSubs/nFolds));
y_act=nan(nIter,nFolds,floor(nSubs/nFolds));

for ii =1:nIter
    fprintf(['Iteration ',num2str(ii),' out of ',num2str(nIter),'\n'])
    
    subIdx=shuffle(subs); 
    subs_temp=subIdx(1:floor(nSubs/nFolds)*nFolds);
    foldIdx=reshape(subs_temp(1:(floor(nSubs/nFolds)*nFolds)),nFolds,floor(nSubs/nFolds));
    
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
                [r,p]=corr(squeeze(fc_both(trainSubs,ed1,ed2)),behav(trainSubs),'type','Spearman');
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
            summary_features(i,1)=nanmean(thisSub(squeeze(importantEdges)>0))-nanmean(thisSub(squeeze(importantEdges)<0));
        end
       
        %train across subjects for each task separately
        mdl=robustfit(summary_features(trainSubs,1),behav(trainSubs));
        
        %test held out subjects
        y_actual=behav(testSub);
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
r_shuff=tanh(z_shuff);

p_para=(1+sum(r_shuff>=mean_r))/(1+nIter);
mse_ave=nanmean(mse_across);
p_para_mse=(1+sum(mse_shuff<=mean(mse_across)))/(1+nIter);

%% plot the results 
lim_bar=600;

figure; set(gcf,'Position',[57 1000 1600 1000]); hold on;
histogram(r_shuff,'FaceColor',[.65 .65 .65],'EdgeAlpha',0)%histogram of shuffled results
histogram(r_across,'FaceColor',[123,50,148]/255,'EdgeAlpha',0,'FaceAlpha',.85)
% histogram(r_shuff,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5],'LineWidth',5)%histogram of shuffled results
bar(mean_r,lim_bar,'FaceColor',[0 0 0],'BarWidth',0.005,'FaceAlpha',1,'EdgeAlpha',0);%bar of the mean r value %[75,0,130]/255
xlim([-0.3,0.4])
ylim([0 lim_bar])
% title({'Chicago-site internal validation',['r=',num2str(median_r),', p=',num2str(p_para)],['mse=',num2str(mse_ave),', p=',num2str(p_para_mse)]})
set(gca,'FontSize',45,'LineWidth',2)
box off
ylabel('Count')
xlabel('r value')
legend('null distribution','observed model performance','Location','northwest');legend('boxoff')


%% Plot a scatter plot from the median r 
%calculate r values for each individual iteration & fold
r_vals=nan(nIter,nFolds);
pred_r_square=nan(nIter,nFolds);
null_pred_r_square=nan(nIter,nFolds);

for ii=1:nIter
    for f=1:nFolds
        r_vals(ii,f)=corr(squeeze(y_pre(ii,f,:)),squeeze(y_act(ii,f,:)),'Type','Pearson');
        
        Y=squeeze(y_act(ii,f,:));
        Yhat=squeeze(y_pre(ii,f,:));
        ssTotal=(Y-nanmean(Y)).^2;
        if sum(ssTotal==0)
            fprintf('Pred: Zerooooooo \n')
        end
        ssError=(Y-Yhat).^2;
        pred_r_square(ii,f)=nanmean(1-ssError./ssTotal);
        
        
        Y_shuff=shuffle(squeeze(y_act(ii,f,:)));
        Yhat=squeeze(y_pre(ii,f,:));
        ssTotal=(Y_shuff-nanmean(Y_shuff)).^2;
        if sum(ssTotal==0)
            fprintf('Null: Zerooooooo \n')
        end
        ssError=(Y_shuff-Yhat).^2;
        null_pred_r_square(ii,f)=nanmean(1-ssError./ssTotal);
    end
end
%calculate the median r value
median_r=nanmedian(r_vals(:));
%find the iteration/fold that is closest to this r value
smallest_dif=abs(r_vals-median_r);
[m,i]=min(smallest_dif);
[mm ii]=min([smallest_dif(i(1),1),smallest_dif(i(2),2),smallest_dif(i(3),3),smallest_dif(i(4),4),smallest_dif(i(5),5)]);

%use the iteration/fold with the median r value as the example for the plot
thisPred=squeeze(y_pre(i(ii),ii,:));
thisActual=squeeze(y_act(i(ii),ii,:));

figure; set(gcf,'Position',[57 1000 1000 800]); hold on
scatter(thisActual,thisPred,200,[0 0 0],'filled','MarkerFaceAlpha',.5)
x=0:5; 
[r_pred,p_pred]=corr(thisPred,thisActual,'Type','Spearman');
[y_pred2,delta_pred]=polyfit(thisActual,thisPred,1);
y_fit=y_pred2(1)*x+y_pred2(2);
plot(x,y_fit,'LineWidth',5,'Color',[0 0 0]);
% title({'Train/test Oregon',['Median r = ',num2str(r_pred)],['MSE = ',num2str(mse)]})
ylabel('Predicted K')
xlabel('Observed K')
set(gca,'FontSize',45,'LineWidth',2)
ylim([0 4]); xlim([0 4])


%% Finished, and let me know
load handel
sound(y,Fs)