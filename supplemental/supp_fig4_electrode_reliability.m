%% Is the reliability of electrodes related to their significance? 
clear all; close all
%% Load the data 

dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/';

%load the data for oregon 
oregon=load([dir,'data/compiled/oregon_raw_singleTrial.mat']);
oregon_elec_idx=load([dir,'data/compiled/oregon_electrode_order.mat']);

load([dir,'/data/compiled/oregon_site_fc_mat.mat'])
%load the edge information
edges_chicago=load([dir,'/data/compiled/edgeWeightsChicago.mat']);
edges_oregon=load([dir,'/data/compiled/edgeWeightsChicago.mat']);

%Load the data for chicago 
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp12_moselle2.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp11_moselle1.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp10_relint2.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp9_relint1.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp8_unpub2.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp7_unpub1.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp6_jocn2.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp5_jocn1.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp4_psychSci4.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp3_psychSci3.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp2_psychSci2.mat'])
load([dir,'/data/compiled/chicago_individual_studiess/amplitude_exp1_psychSci1.mat'])

nIter=10000;
nElec=length(oregon_elec_idx.bs_e_idx); 

%% Oregon: Split half reliability of each electrode across trials for each subject
tmp_oregon_trial=squeeze(nanmean(oregon.data.contraIpsiCorr.allTrial(:,:,oregon_elec_idx.bs_e_idx,ismember(-200:4:999,0:1000)),4));

oregon.split_reliability.trials=nan(nIter,nElec);
for ii=1:nIter
    fprintf(['Oregon iteration ',num2str(ii),' out of ',num2str(nIter),'\n'])
    for e=1:nElec
        trialIdx=shuffle(1:size(tmp_oregon_trial,2));
        nHalf=round(length(trialIdx)/2);
        oregon.split_reliability.trials(ii,e)=corr(squeeze(nanmean(tmp_oregon_trial(:,trialIdx(1:nHalf),e),2)),squeeze(nanmean(tmp_oregon_trial(:,trialIdx(nHalf+1:end),e),2)));
    end
end

%% Chicago: Cronbach's alpha per electrode, across trials 
tmp_chicago_trials=cat(1,psychSci1.amplitude_all_trials,psychSci2.amplitude_all_trials,psychSci3.amplitude_all_trials,...
    psychSci4.amplitude_all_trials,jocn1.amplitude_all_trials,jocn2.amplitude_all_trials,...
    moselle1.amplitude_all_trials,moselle2.amplitude_all_trials,relint1.amplitude_all_trials,...
    relint2.amplitude_all_trials,unpub1.amplitude_all_trials,unpub2.amplitude_all_trials);

chicago.split_reliability.trials=nan(nIter,nElec); 
for ii=1:nIter
    fprintf(['Chicago iteration ',num2str(ii),' out of ',num2str(nIter),'\n'])
    for e=1:nElec
        trialIdx=shuffle(1:size(tmp_chicago_trials,2));
        nHalf=round(length(trialIdx)/2);
        
        chi_1=squeeze(nanmean(tmp_chicago_trials(:,trialIdx(1:nHalf),e),2));
        chi_2=squeeze(nanmean(tmp_chicago_trials(:,trialIdx(nHalf+1:end),e),2));
        
        chicago.split_reliability.trials(ii,e)=corr(chi_1,chi_2);
    end
end

%% Does the reliability of each electrode relate to whether it is predictive of behavior? 
[h p]=corr(mean(edges_oregon.r_chi,2),squeeze(nanmean(oregon.split_reliability.trials,1))')

[h p]=corr(nanmean(edges_chicago.r_chi,2),squeeze(nanmean(chicago.split_reliability.trials,1))')

%% ave split-half reliability
oregon_reliability=squeeze(nanmean(oregon.split_reliability.trials,1));
chicago_reliabillity=squeeze(nanmean(chicago.split_reliability.trials,1));

%% Plot split-half reliability results (for Supplemental Materials)

figure; 

subplot(2,1,1)
imagesc(oregon_reliability)
h=colorbar; caxis([0.7 1]); ylabel(h, 'r value');
set(gca,'FontSize',30,'LineWidth',2);box off;
set(gca,'XTick',[1:17],'XTickLabel',bs.bs.fc_mat.overlappingElecs.elecs,'TickDir','out')
set(gca,'YTick',[],'YTickLabel',{})
xtickangle(45)
title({'Electrode split-half reliability','Oregon-site data'})

subplot(2,1,2)
imagesc(chicago_reliabillity)
colorbar; caxis([0.7 1])
set(gca,'FontSize',30,'LineWidth',2);box off;
h=colorbar; caxis([0.7 1]); ylabel(h, 'r value');
set(gca,'XTick',[1:17],'XTickLabel',bs.bs.fc_mat.overlappingElecs.elecs,'TickDir','out')
set(gca,'YTick',[],'YTickLabel',{})
xtickangle(45)
xlabel('electrode')
title('Chicago-site data')

%% Figure v2

figure 
imagesc([oregon_reliability;chicago_reliabillity])
h=colorbar; ylabel(h, 'r value');
set(gca,'FontSize',30,'LineWidth',2);box off;
set(gca,'XTick',[1:17],'XTickLabel',bs.bs.fc_mat.overlappingElecs.elecs,'TickDir','out')
set(gca,'YTick',[1,2],'YTickLabel',{'Oregon','Chicago'})
xtickangle(45)
xlabel('electrode')
title('Electrode reliability')

%% Figure v3 (lines)

figure; hold on 
plot(oregon_reliability,'linewidth',6,'color',[59 46 102]/255);
plot(chicago_reliabillity,'linewidth',6,'color',[135 105 234]/255);
set(gca,'XTick',[1:17],'XTickLabel',bs.bs.fc_mat.overlappingElecs.elecs,'TickDir','out')
xtickangle(45)
xlabel('electrode')
title('Electrode reliability')
set(gca,'FontSize',30,'LineWidth',2);box off;
legend('Oregon','Chicago')



