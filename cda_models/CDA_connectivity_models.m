dbstop if error;
clear all
%% Load the Oregon data 
%raw, filtered EEG data
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/BigStudy/Data/BigStudy_filtered_erp_singleTrial.mat')
%fc_mat data
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/BigStudy/Data/fcmat_allsubs_RI_and_baseline_filtered.mat')
%load important electrodes
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/oregon_chicago_top10percedges.mat')
%% Must get rid of extra subbies
times=-200:4:999; 
cda_times=300:999;
cda_time_idx=ismember(times,cda_times);
nIter=10000;

eeg_idx=ismember(data.subs,beh.cd_cda.subNum);
importantElectrodes={'F3'    'F4'    'C3'    'C4'    'P7'    'P8'    'P3'    'P4'    'PO7'    'PO8'    'PO3'    'PO4'    'O1'    'O2'    'Fz'    'Cz'    'Pz'};
elec_idx=nan(length(importantElectrodes),1);
for elec=1:length(importantElectrodes)
   elec_idx(elec)=find(contains(data.tempEEGelec,importantElectrodes(elec)));
end

%contra/ipsi
contra_ipsi_elecs={'PO3','PO4';'F3','F4';'C3','C4';'P3','P4';'O1','O2';'PO7','PO8';'P7','P8'};
contra_ipsi_idx=nan(size(contra_ipsi_elecs,1),size(contra_ipsi_elecs,2));
for elec=1:size(contra_ipsi_elecs,1)
    for elec2=1:size(contra_ipsi_elecs,2)
        contra_ipsi_idx(elec,elec2)=find(contains(data.tempEEGelec,contra_ipsi_elecs(elec,elec2)));
    end
end

oregon.amplitude_data.color=squeeze(nanmean(cat(2,data.contraIpsiCorr.C2(eeg_idx,:,elec_idx,:),data.contraIpsiCorr.C6(eeg_idx,:,elec_idx,:)),2));%get all electrodes
oregon.amplitude_data.color=oregon.amplitude_data.color(:,contra_ipsi_idx(:,1),:)-oregon.amplitude_data.color(:,contra_ipsi_idx(:,2),:);%then subtract contra-ipsi for the cda electrodes
oregon.cda.color=squeeze(nanmean(nanmean(oregon.amplitude_data.color(:,[1,4:7],cda_time_idx),2),3));%then, only get the relevant times and electrodes and then average, so that there is only one value per person

oregon.amplitude_data.shape=squeeze(nanmean(cat(2,data.contraIpsiCorr.S2(eeg_idx,:,elec_idx,:),data.contraIpsiCorr.S6(eeg_idx,:,elec_idx,:)),2));
oregon.amplitude_data.shape=oregon.amplitude_data.shape(:,contra_ipsi_idx(:,1),:)-oregon.amplitude_data.shape(:,contra_ipsi_idx(:,2),:);%then subtract contra-ipsi for the cda electrodes
oregon.cda.shape=squeeze(nanmean(nanmean(oregon.amplitude_data.shape(:,[1,4:7],cda_time_idx),2),3));%then, only get the relevant times and electrodes and then average, so that there is only one value per person

oregon.cda.both=nanmean([oregon.cda.shape,oregon.cda.color],2);%this is what we will use for the CDA for Oregon

oregon.fc.color=bs.fc_mat.overlappingElecs.RI.color; 
oregon.fc.shape=bs.fc_mat.overlappingElecs.RI.shape; 

nSubs.oregon=size(oregon.fc.color,1);

%% Load the Chicago data 
chicago=load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/fc_mat_uchicagosubs_n165.mat');%'/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/fc_mat_matchedTrials_uchicagosubs_n165_vJune2020.mat')

chicago.cda=nanmean(chicago.fc_mat.contraMipsi(:,[1,4:7]),2);%get the cda electrodes and then average

nSubs.chicago=length(chicago.cda);
%% Get the summary feature for each subject - Oregon
oregon.summary_features.color=nan(nSubs.oregon,1);
oregon.summary_features.shape=nan(nSubs.oregon,1);
chicago.summary_features.chicago_color=nan(nSubs.oregon,1);
chicago.summary_features.chicago_shape=nan(nSubs.oregon,1);
for s=1:nSubs.oregon
    thisSub_shape=squeeze(oregon.fc.shape(s,:,:));
    thisSub_color=squeeze(oregon.fc.color(s,:,:));
    
    oregon.summary_features.shape(s)=nanmean(thisSub_shape(sigEdges.shape==1))-nanmean(thisSub_shape(sigEdges.shape==-1));
    oregon.summary_features.color(s)=nanmean(thisSub_color(sigEdges.color==1))-nanmean(thisSub_color(sigEdges.color==-1));
    
    chicago.summary_features.chicago_color(s)=nanmean(thisSub_color(sigEdges.chicago==1))-nanmean(thisSub_color(sigEdges.chicago==-1));
    chicago.summary_features.chicago_shape(s)=nanmean(thisSub_shape(sigEdges.chicago==1))-nanmean(thisSub_shape(sigEdges.chicago==-1));
end


oregon.summary_features.trainOregon=nanmean([oregon.summary_features.shape,oregon.summary_features.color],2);%this is what we will use for connectivity for Oregon
oregon.summary_features.trainChicago=nanmean([chicago.summary_features.chicago_color,chicago.summary_features.chicago_shape],2);
%% Get the summary features for each subject - Chicago
chicago.summary_features.trainChicago=nan(nSubs.chicago,1); 
chicago.summary_features.trainOregon=nan(nSubs.chicago,1);

for s=1:nSubs.chicago
   thisSub=squeeze(chicago.fc_mat.eeg(s,:,:));
   chicago.summary_features.trainChicago(s)=nanmean(thisSub(sigEdges.chicago==1))-nanmean(thisSub(sigEdges.chicago==-1));
   chicago.summary_features.trainOregon(s)=nanmean(thisSub(sigEdges_overlap.color_shape==1))-nanmean(thisSub(sigEdges_overlap.color_shape==-1));
end
%% CDA models

%train the models
mdl_cda.trainOregon=fitlm(oregon.cda.both,bs.fc_mat.k.ave);
mdl_cda.trainChicago=fitlm(chicago.cda,chicago.fc_mat.beh);

%make predictions
pred_cda.trainOregon=predict(mdl_cda.trainOregon,chicago.cda);
pred_cda.trainChicago=predict(mdl_cda.trainChicago,oregon.cda.both);

%calculate shuffled r values 
r_trainO_cda=corr(pred_cda.trainOregon,chicago.fc_mat.beh,'type','Spearman');
r_trainC_cda=corr(pred_cda.trainChicago,bs.fc_mat.k.ave,'type','Spearman');
r_shuff_trainO_cda=nan(1,nIter);
r_shuff_trainC_cda=nan(1,nIter);
for ii=1:nIter
    r_shuff_trainO_cda(ii)=corr(shuffle(pred_cda.trainOregon),chicago.fc_mat.beh,'type','Spearman');
    r_shuff_trainC_cda(ii)=corr(shuffle(pred_cda.trainChicago),bs.fc_mat.k.ave,'type','Spearman');
end
p_trainO_cda=(1+sum(r_shuff_trainO_cda>=r_trainO_cda))/(1+nIter);
p_trainC_cda=(1+sum(r_shuff_trainC_cda>=r_trainC_cda))/(1+nIter);

%calculate mse & shuffled mse & significance
mse_trainO_cda=immse(pred_cda.trainOregon,chicago.fc_mat.beh);
mse_trainC_cda=immse(bs.fc_mat.k.ave,pred_cda.trainChicago);
mse_shuff_trainO_cda=nan(1,nIter);
mse_shuff_trainC_cda=nan(1,nIter);
for ii=1:nIter
    mse_shuff_trainO_cda(ii)=immse(shuffle(pred_cda.trainOregon),chicago.fc_mat.beh);
    mse_shuff_trainC_cda(ii)=immse(shuffle(bs.fc_mat.k.ave),pred_cda.trainChicago);
end
p_mse_traino_cda=(1+sum(mse_shuff_trainO_cda<=mse_trainO_cda))/(1+nIter);
p_mse_trainc_cda=(1+sum(mse_shuff_trainC_cda<=mse_trainC_cda))/(1+nIter);

%plot the results of these two models 
figure; set(gcf,'Position',[57 1000 2300 1000]);

%train oregon, test chicago
subplot(1,2,1)
scatter(chicago.fc_mat.beh,pred_cda.trainOregon,200,'filled','MarkerFaceAlpha',.4)
[h,p]=corr(pred_cda.trainOregon,chicago.fc_mat.beh,'type','Spearman');
title({'Train Oregon, test Chicago','','CDA model',['r=',num2str(h),', p=',num2str(p)],['mse=',num2str(mse_trainO_cda),', p=',num2str(p_mse_traino_cda)]})
xlabel('Actual K')
ylabel('Predicted K')
set(gca,'FontSize',20,'LineWidth',3)
xlim([0 5]);ylim([0 6])

%train chicago, test oregon
subplot(1,2,2)
scatter(bs.fc_mat.k.ave,pred_cda.trainChicago,200,'filled','MarkerFaceAlpha',.4);
[h,p]=corr(pred_cda.trainChicago,bs.fc_mat.k.ave,'type','Spearman');
title({'Train chicago, test Oregon','','CDA model',['r=',num2str(h),', p=',num2str(p)],['mse=',num2str(mse_trainC_cda),', p=',num2str(p_mse_trainc_cda)]})
xlabel('Actual K')
ylabel('Predicted K')
set(gca,'FontSize',20,'LineWidth',3)
xlim([0 5]);ylim([0 6])

figure; 
scatter(bs.fc_mat.k.ave,pred_cda.trainChicago,200,'filled','MarkerFaceAlpha',.4);
[h,p]=corr(pred_cda.trainChicago,bs.fc_mat.k.ave,'type','Spearman');
title({'Train chicago, test Oregon','','CDA model',['r=',num2str(h),', p=',num2str(p),', mse=',num2str(immse(bs.fc_mat.k.ave,pred_cda.trainChicago))]})
xlabel('Actual K')
ylabel('Predicted K')
set(gca,'FontSize',20,'LineWidth',3)
xlim([0 5]);ylim([2.46 2.54])

%% FC models 

%train the models
mdl_fc.trainOregon=fitlm(oregon.summary_features.trainOregon,bs.fc_mat.k.ave);
mdl_fc.trainChicago=fitlm(chicago.summary_features.trainChicago,chicago.fc_mat.beh);

%make predictions
pred_fc.trainOregon=predict(mdl_fc.trainOregon,chicago.summary_features.trainOregon);
pred_fc.trainChicago=predict(mdl_fc.trainChicago,oregon.summary_features.trainChicago); 

%train oregon, test chicago
subplot(3,2,3)
scatter(chicago.fc_mat.beh,pred_fc.trainOregon,200,'filled','MarkerFaceAlpha',.4)
[h,p]=corr(chicago.fc_mat.beh,pred_fc.trainOregon,'type','Spearman');
title({'FC model',['r=',num2str(h),', p=',num2str(p)]})
xlabel('Actual K')
ylabel('Predicted K')
set(gca,'FontSize',20,'LineWidth',3)
xlim([0 5]);ylim([0 3.2])

%train chicago, test oregon
subplot(3,2,4)
scatter(bs.fc_mat.k.ave,pred_fc.trainChicago,200,'filled','MarkerFaceAlpha',.4)
[h,p]=corr(bs.fc_mat.k.ave,pred_fc.trainChicago,'type','Spearman');
title({'FC model',['r=',num2str(h),', p=',num2str(p)]})
xlabel('Actual K')
ylabel('Predicted K')
set(gca,'FontSize',20,'LineWidth',3)
xlim([0 5]);ylim([0 3.2])

%% FC & CDA model (combo)

%train the models 
mdl_both.trainOregon=fitlm([oregon.summary_features.trainOregon,oregon.cda.both],bs.fc_mat.k.ave);
mdl_both.trainChicago=fitlm([chicago.summary_features.trainChicago,chicago.cda],chicago.fc_mat.beh);

%make predictions
pred_both.trainOregon=predict(mdl_both.trainOregon,[chicago.summary_features.trainOregon,chicago.cda]);
pred_both.trainChicago=predict(mdl_both.trainChicago,[oregon.summary_features.trainChicago,oregon.cda.both]);


%calculate shuffled r values 
r_trainO_both=corr(pred_both.trainOregon,chicago.fc_mat.beh,'type','Spearman');
r_trainC_both=corr(pred_both.trainChicago,bs.fc_mat.k.ave,'type','Spearman');
r_shuff_trainO_both=nan(1,nIter);
r_shuff_trainC_both=nan(1,nIter);
for ii=1:nIter
    r_shuff_trainO_both(ii)=corr(shuffle(pred_both.trainOregon),chicago.fc_mat.beh,'type','Spearman');
    r_shuff_trainC_both(ii)=corr(shuffle(pred_both.trainChicago),bs.fc_mat.k.ave,'type','Spearman');
end
p_trainO_both=(1+sum(r_shuff_trainO_both>=r_trainO_both))/(1+nIter);
p_trainC_both=(1+sum(r_shuff_trainC_both>=r_trainC_both))/(1+nIter);

%calculate mse & shuffled mse & significance
mse_trainO_both=immse(pred_both.trainOregon,chicago.fc_mat.beh);
mse_trainC_both=immse(bs.fc_mat.k.ave,pred_both.trainChicago);
mse_shuff_trainO_both=nan(1,nIter);
mse_shuff_trainC_both=nan(1,nIter);
for ii=1:nIter
    mse_shuff_trainO_both(ii)=immse(shuffle(pred_both.trainOregon),chicago.fc_mat.beh);
    mse_shuff_trainC_both(ii)=immse(shuffle(bs.fc_mat.k.ave),pred_both.trainChicago);
end
p_mse_traino_both=(1+sum(mse_shuff_trainO_both<=mse_trainO_both))/(1+nIter);
p_mse_trainc_both=(1+sum(mse_shuff_trainC_both<=mse_trainC_both))/(1+nIter);


%train oregon, test chicago
subplot(3,2,5)
scatter(chicago.fc_mat.beh,pred_both.trainOregon,200,'filled','MarkerFaceAlpha',.4)
[h,p]=corr(chicago.fc_mat.beh,pred_both.trainOregon,'type','Spearman');
title({'FC + CDA model',['r=',num2str(h),', p=',num2str(p)],['mse=',num2str(mse_traino_combo),', p=',num2str(p_mse_traino)]})
xlabel('Actual K')
ylabel('Predicted K')
set(gca,'FontSize',20,'LineWidth',3)
xlim([0 5]);ylim([0 4.5])

%train chicago, test oregon 
subplot(3,2,6)
scatter(bs.fc_mat.k.ave,pred_both.trainChicago,200,'filled','MarkerFaceAlpha',.4)
[h,p]=corr(bs.fc_mat.k.ave,pred_both.trainChicago,'type','Spearman');
title({'FC + CDA model',['r=',num2str(h),', p=',num2str(p)],['mse=',num2str(mse_trainc_combo),', p=',num2str(p_mse_trainc)]})
xlabel('Actual K')
ylabel('Predicted K')
set(gca,'FontSize',20,'LineWidth',3)
xlim([0 5]);ylim([0 4.5])

%% Plot correlation between CDA during CDA time window

figure; set(gcf,'Position',[57 1000 2300 1000]);

subplot(1,2,1)
scatter(oregon.cda.both,bs.fc_mat.k.ave,200,'filled','MarkerFaceAlpha',.4)
[h,p]=corr(oregon.cda.both,bs.fc_mat.k.ave,'type','Spearman');
title({'Correlation between CDA and K score','Oregon',['r=',num2str(h),', p=',num2str(p)]})
xlabel('CDA');ylabel('K score')
set(gca,'FontSize',20,'LineWidth',3)
xlim([-9 3]);ylim([0 5])

subplot(1,2,2)
scatter(chicago.cda,chicago.fc_mat.beh,200,'filled','MarkerFaceAlpha',.4)
[h,p]=corr(chicago.cda,chicago.fc_mat.beh,'type','Spearman');
title({'Chicago',['r=',num2str(h),', p=',num2str(p)]})
xlabel('CDA');ylabel('K score')
set(gca,'FontSize',20,'LineWidth',3)
xlim([-9 3]);ylim([0 5])

%% Plot correlation between CDA and beahvior (with each chicago exp colored differently) 
uniqColors=distinguishable_colors(12);% unique colors 

figure; 

subplot(2,1,1)
for e=1:12
    hold on
    exp_ind=logical(sum(ismember(chicago.fc_mat.experiment,e),2)); 
%     scatter(chicago.cda(exp_ind),chicago.fc_mat.beh(exp_ind),100,uniqColors(e,:),'filled','MarkerFaceAlpha',.4)
    
    x=-10:5;
    [y_pred,delta_pred]=polyfit(chicago.cda(exp_ind),chicago.fc_mat.beh(exp_ind),1);
    y_fit=y_pred(1)*x+y_pred(2);
    plot(x,y_fit,'LineWidth',5,'Color',uniqColors(e,:));
end
[h,p]=corr(chicago.cda,chicago.fc_mat.beh,'type','Spearman');
legend('PsycSci1','PsycSci2','PsycSci3','PsycSci4',...
    'JoCN1','JoCN2','UnpubInt1','UnpubInt2','relInt1','relInt2','moselle1','moselle2')
xlim([-9 3]);ylim([0 5])
xlabel('CDA');ylabel('K score')
set(gca,'FontSize',20,'LineWidth',3)
title({'Correlation between CDA and K score','Chicago: separated by experiment',['r=',num2str(h),', p=',num2str(p)]})

subplot(2,1,2)
for e=1:12
    hold on
    exp_ind=logical(sum(ismember(chicago.fc_mat.experiment,e),2)); 
    scatter(chicago.cda(exp_ind),chicago.fc_mat.beh(exp_ind),100,uniqColors(e,:),'filled','MarkerFaceAlpha',.4)
    
end
xlim([-9 3]);ylim([0 5])
xlabel('CDA');ylabel('K score')
set(gca,'FontSize',20,'LineWidth',3)

%% Scatter plot of CDA for both Oregon and Chicago 
figure;set(gcf,'Position',[57 1000 1000 700]);
beh=[oregon.cda.both,[chicago.cda;nan(6,1)]];
beh=num2cell(beh,1); 
 
BF_JitteredParallelScatter_MdB_wm(beh,1,1,0)
errorbar([1:2],[nanmean(oregon.cda.both),nanmean(chicago.cda)],[nanstd(oregon.cda.both)/sqrt(length(oregon.cda.both)),nanstd(chicago.cda)/sqrt(length(chicago.cda))],'.','LineWidth',8,'Color',[0 0 0],'CapSize',0)%error bars

set(gca,'XTick',[1:1:2],'XTickLabel',{['Oregon: m=',num2str(nanmean(oregon.cda.both)),' sem=',num2str(nanstd(oregon.cda.both)/sqrt(length(oregon.cda.both)))],...
    ['Chicago: m=',num2str(nanmean(chicago.cda)),' sem=',num2str(nanstd(chicago.cda)/sqrt(length(chicago.cda)))]})
set(gca,'FontSize',25,'LineWidth',2,'box','off','TickDir','out')
ylabel('CDA amplitude')
title('CDA amplitude for Oregon vs. Chicago')

