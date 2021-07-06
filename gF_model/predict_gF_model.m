%% Train and test gF with 5-fold CV, using its own edges
clear all; close all; 
%% Load the data 
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/BigStudy/Data/fcmat_allsubs_RI_and_baseline_filtered.mat')%Oregon
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/BigStudy/Data/all_tasks_dat.mat')
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/BigStudy/Data/BigStudy_filtered_erp_singleTrial.mat')%need to load this for the EEG subject numbers

% do these analyses with the edges that overlapped 
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/oregon_chicago_top10percedges.mat')%sig_edges_trainOreg_testChi.mat')
importantEdges=sigEdges_overlap.oregon_chicago;%do the  only with the overlapping edges from Oregon and Chicago for both datasets

eeg_subs=data.subs(ismember(data.subs,beh.cd_cda.subNum));
eeg_sub_idx=ismember(eeg_subs,dat.subjects); 
behav_other_subs=ismember(dat.subjects,eeg_subs);

fc_both    = mean(cat(4,bs.fc_mat.overlappingElecs.RI.color(eeg_sub_idx,:,:),bs.fc_mat.overlappingElecs.RI.shape(eeg_sub_idx,:,:)),4);% big_study_fc    = mean(cat(4,bs.fc_mat.color,bs.fc_mat.shape),4);
behav = dat.gF(behav_other_subs);

nElec  = size(fc_both,2);%17;   % number of nodes

nIter=10000; 
nFolds=5; 
subs=1:sum(behav_other_subs);
nSubs=length(subs); 

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
median_r=tanh(median(z_across));
r_shuff=tanh(z_shuff);

p_para=(1+sum(r_shuff>=mean_r))/(1+nIter);
mse_ave=nanmean(mse_across);
p_para_mse=(1+sum(mse_shuff<=mean(mse_across)))/(1+nIter);

%% Plot a scatter plot from the median r 
%% plot the results 
lim_bar=900;

figure; set(gcf,'Position',[57 1000 1600 1000]); hold on;
histogram(r_shuff,'FaceColor',[.65 .65 .65],'EdgeAlpha',0)%histogram of shuffled results
histogram(r_across,'FaceColor',[123,50,148]/255,'EdgeAlpha',0,'FaceAlpha',.85)
histogram(r_shuff,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5],'LineWidth',5)%histogram of shuffled results
bar(median_r,lim_bar,'FaceColor',[0 0 0],'BarWidth',0.005,'FaceAlpha',1,'EdgeAlpha',0);%nanmean(r_across),lim_bar,'FaceColor',[0 0 0],'BarWidth',0.005,'FaceAlpha',1,'EdgeAlpha',0)%bar of the mean r value %[75,0,130]/255
xlim([-0.3,0.4])
ylim([0 lim_bar])
set(gca,'FontSize',45,'LineWidth',2)
box off
ylabel('Count')
xlabel('r value')
legend('null distribution','observed model performance','Location','northwest');legend('boxoff')
%% Plot median r

%calculate r values for each individual iteration & fold
r_vals=nan(nIter,nFolds);
for ii=1:nIter
    for f=1:nFolds
        r_vals(ii,f)=corr(squeeze(y_pre(ii,f,:)),squeeze(y_act(ii,f,:)),'Type','Spearman');
    end
end
%calculate the median r value
mean_r=nanmean(r_vals(:));
%find the iteration/fold that is closest to this r value
smallest_dif=abs(r_vals-mean_r);
[m,i]=min(smallest_dif);
[mm ii]=min([smallest_dif(i(1),1),smallest_dif(i(2),2),smallest_dif(i(3),3),smallest_dif(i(4),4),smallest_dif(i(5),5)]);

%use the iteration/fold with the median r value as the example for the plot
thisPred=squeeze(y_pre(i(ii),ii,:));
thisActual=squeeze(y_act(i(ii),ii,:));

x=8:25; 
[r_pred,p_pred]=corr(thisPred,thisActual,'Type','Spearman');
[y_pred2,delta_pred]=polyfit(thisActual,thisPred,1);
y_fit=y_pred2(1)*x+y_pred2(2);
mse=mean((thisActual-thisPred).^2);
figure; set(gcf,'Position',[57 1000 1000 800]); hold on
scatter(thisActual,thisPred,200,[0 0 0],'filled','MarkerFaceAlpha',.5)
plot(x,y_fit,'LineWidth',5,'Color',[0 0 0]);
% title({'Train/test Oregon',['Median r = ',num2str(r_pred)],['MSE = ',num2str(mse)]})
ylabel('Predicted K')
xlabel('Observed K')
set(gca,'FontSize',45,'LineWidth',2)
ylim([12 22]); xlim([12 22])


%% Finished, and let me know
load handel
sound(y,Fs)