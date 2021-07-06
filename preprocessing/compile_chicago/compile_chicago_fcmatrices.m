%% This script imports the subject numbers and unique IDs from the UChicago latearlized change detection studies that I have identified
% This includes data from: 
% 1. Psychological science studies (4 experiments)
% 2. Journal of Cognitive Neuroscience studies (2 experiments) 
% 3. Unpublished interruption experiments, which were part of JoCN data collection (3 experiments)
% 4. Relevant interruption experiments (2 of which will be published, 1 of which will not be published)
% 5. Moselle experiment - run by Tobias Feldmann-Wüstefeld
% February 2020, Nicole Hakim 
%% clear and change directory
clear all; clc
cd('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI')
savePath='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/';
%% Settings to change - which time period do you want to use 
whichTimes=2; %1=baseline, 2=retention interval

if whichTimes==1%baseline time window 
    saveName='fc_mat_baselineTimes_uchicagosubs_n';
elseif whichTimes==2%retention interval time window
    saveName='chicago_site_fc_mat';
end

matchTrials=1;%1 =if we want to match the number of trials across all participants. The participant with the least number of trials had 217 trials, so we will make everyone have that number of trials. We will only include trials from the first experiment that they did

if matchTrials==1
   saveName='fc_mat_matchedTrials_uchicagosubs_n'; 
end
%% Path settings
path_psycSci='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_PsychSci/data/';
path_jocn='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_JoCN/data/';
path_unpubint='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_unpublished_int/data/';
path_relint='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_RelInt/data/';
path_moselle='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Toby_Moselle/data/';
%% Load the subject numbers and unique ids for each experiment 
%psych science
load([path_psycSci,'PsycSci_uniqueIds.mat']);
uniqueSubs.psycSci=unique([psycSci.uniqueIds.e1,psycSci.uniqueIds.e2,psycSci.uniqueIds.e3,psycSci.uniqueIds.e4]');

%jocn
load([path_jocn,'jocn_uniqIds_subNums.mat']);
jocn.subs.uniqIds=[jocn.uniqueIds.e2;jocn.uniqueIds.e1];
uniqueSubs.jocn=unique([jocn.uniqueIds.e2,jocn.uniqueIds.e1]);

%unpublished int
load([path_unpubint,'unpub_int_uniqIds_subNums.mat'])
unpub_int.uniqIds=[unpub_int.uniqueIds.e4,unpub_int.uniqueIds.e5];
uniqueSubs.unpub=unique([unpub_int.uniqueIds.e4,unpub_int.uniqueIds.e5]);

%relevant int
load([path_relint,'relint_uniqIds_subNums.mat']); 
uniqueSubs.relint=unique([relint.uniqueIds.e1,relint.uniqueIds.e2]);

%moselle 
load([path_moselle,'moselle_uniqIds_subNums.mat'])
uniqueSubs.moselle=unique([moselle.uniqueIds.e1,moselle.uniqueIds.e2]);

%figure out the unique subjects in all experiments 
uniqueSubs.all=unique([uniqueSubs.psycSci',uniqueSubs.jocn,uniqueSubs.unpub,uniqueSubs.relint,uniqueSubs.moselle]);
nUniqueSubs=length(uniqueSubs.all);

%% Get the times for each experiment -- this is for the retention interval

%save the time range to a structure for each experiment
if whichTimes==2
    Oregon_RI_length=1000;
    maxTimes=min([Oregon_RI_length,psycSci.times.e1(end),jocn.times.e2(end),jocn.times.e1(end),unpub_int.times.e5(end),relint.times.e2(end),moselle.times.e2(end)]);
    timeIdx=0:4:maxTimes;%by 4 because the Oregon data is sampled at 250 Hz, which is the lowest of all of the datasets 
end
%% If we are interested in the baseline time window, then load some data about the basline lengths of each experiment and then find which experiment had the shortest baseline window. 
if whichTimes==1
    load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_JoCN/data/experiment2/6_eeg_beh_jocn_e2.mat')
    jocn_baslineTimes.e2=jocn_e2.times_baseline;
    
    load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_unpublished_int/data/experiment4/9_eeg_beh_unpubint_e4.mat')
    unpubint_baselineTimes.e4=unpubint_e4.times_baseline;
    
    load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_RelInt/data/experiment1/2_eeg_beh_relint_e1.mat')
    relint_baselineTimes.e1=relint_e1.times_baseline;
    
    load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_RelInt/data/experiment2/2_eeg_beh_relint_e2.mat')
    relint_baselineTimes.e2=relint_e2.times_baseline;
    
    load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Toby_Moselle/data/experiment1/13_eeg_beh_moselle_e1.mat')
    moselle_baselineTimes.e1=moselle_e1.times_baseline;
    
    % Get the times for each experiment -- this is for the baseline period
    Oregon_pre_trial_baseline=-50;
    maxTimes_base=max([Oregon_pre_trial_baseline,psycSci.times_baseline.e1(1),jocn_baslineTimes.e2(1),...
        unpubint_baselineTimes.e4(1),relint_baselineTimes.e1(1),relint_baselineTimes.e2(1),moselle_baselineTimes.e1(1)]);
    timeIdx=maxTimes_base:4:0;%by 4 because the Oregon data is sampled at 250 Hz, which is the lowest of all of the datasets 
end
%% Order of chicago electrodes
raw_order={'F3','F4','C3','C4','P7','P8','P3','P4',...
    'PO7','PO8','PO3','PO4','O1','O2','Fz','Cz','Pz'};
% save the raw contra - ipsi electrodes
contra_ipsi_elecs={'PO3','PO4';'F3','F4';'C3','C4';'P3','P4';'O1','O2';'PO7','PO8';'P7','P8'};
%% 
psycSci.subs.nums(4,27)=37;% this will be fixed in the future
%% Loop through subjects & save their data into one large fc matrix 
nElec=17; %number of overlapping electrodes between uchicago and oregon

fc_mat.eeg=nan(nUniqueSubs,nElec,nElec);                                 %functional connectivity matrix
fc_mat.eeg_allE=nan(nUniqueSubs,nElec,nElec);                            %functional connectivity matrix - all 30 electrodes 
fc_mat.contraMipsi=nan(nUniqueSubs,size(contra_ipsi_elecs,1));           %contra minus ipsi for posterior electrodes
contraMipsi_400_to_1000=nan(nUniqueSubs,size(contra_ipsi_elecs,1));      %contra minus ipsi for posterior electrodes (from 400-1000 ms)
fc_mat.beh=nan(nUniqueSubs,1);                                           %behavioral measures 
fc_mat.uniqueIds=uniqueSubs.all;                                         %unique ids
fc_mat.experiment=nan(nUniqueSubs,12);
fc_mat.elecOrder=raw_order; 
fc_mat.contraIpsiOrder=contra_ipsi_elecs;
fc_mat.nTrials=nan(nUniqueSubs,1);  
fc_mat.amplitude.raw=nan(nUniqueSubs,nElec,251);
fc_mat.amplitude.raw_aveTime=nan(nUniqueSubs,nElec);

for s=1:nUniqueSubs
    tic
    fprintf(['Subject ',num2str(s),' out of ',num2str(nUniqueSubs),'\n'])
    
    thisId=uniqueSubs.all(s);%this participant's unique ID
    
    %which experiments did this participant do? 
    whichExp=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
        ismember(thisId,psycSci.uniqueIds.e2),...        %Psych sci exp 2
        ismember(thisId,psycSci.uniqueIds.e3),...        %Psych sci exp 3
        ismember(thisId,psycSci.uniqueIds.e4),...        %Psych sci exp 4
        ismember(thisId,jocn.uniqueIds.e1),...                  %JoCN exp 1
        ismember(thisId,jocn.uniqueIds.e2),...                  %JoCN exp 2
        ismember(thisId,unpub_int.uniqueIds.e4),...             %unpublished int exp 4
        ismember(thisId,unpub_int.uniqueIds.e5),...             %unpublished int exp 5
        ismember(thisId,relint.uniqueIds.e1),...                %relevant interruption experiment 1
        ismember(thisId,relint.uniqueIds.e2),...                %relevant interruption experiment 2
        ismember(thisId,moselle.uniqueIds.e1),...               %moselle experiment 1
        ismember(thisId,moselle.uniqueIds.e2)];                 %moselle experiment 2
    
    tmpeeg={};%temporary structure to hold the relevant eeg data 
    tmpeeg_all={};%temorary structure to hold the relevant eeg data (all 30 electrodes)
    tmpbeh=nan(1,12);%temporary structure to hold the relevant behavioral data 
    tmpbeh_wo_ss2=nan(1,12);%temporary structure to hold k score (not including set size 2)
    tmpbeh_only_ss2=nan(1,12);%temporary structure to hold k score (only including set size 2)
    %loop through the experiments & grab the relevant data for this participant
    if whichExp(1)==1%if this participant did psych sci exp 1
        tmpnum=psycSci.sub_nums.e1(psycSci.uniqueIds.e1==thisId);
        load([path_psycSci,'/experiment1/',num2str(tmpnum),'_eeg_beh_psycsci_e1.mat'])
        if whichTimes==1
            tmpeeg{1}=psychsci_e1.raw_baseline(:,:,ismember(psycSci.times_baseline.e1,timeIdx));
        elseif whichTimes==2
            tmpeeg{1}=psychsci_e1.raw(:,:,ismember(psychsci_e1.times,timeIdx));%get only the times that align for all experiments
            tmpeeg_all{1}=psychsci_e1.raw_allE(:,:,ismember(psychsci_e1.times,timeIdx));%get only the times that align for all experiments
        end
        tmpbeh(1)=psychsci_e1.k;
        tmpbeh_wo_ss2(1)=psychsci_e1.k_ss4only;
        tmpbeh_only_ss2(1)=psychsci_e1.k_ss2only;
    end
    if whichExp(2)==1%if this participant did psych sci exp 2
        tmpnum=psycSci.sub_nums.e2(psycSci.uniqueIds.e2==thisId);
        load([path_psycSci,'/experiment2/',num2str(tmpnum),'_eeg_beh_psycsci_e2.mat'])
        if whichTimes==1
            tmpeeg{2}=psychsci_e2.raw_baseline(:,:,ismember(psycSci.times_baseline.e2,timeIdx));
        elseif whichTimes==2
            tmpeeg{2}=psychsci_e2.raw(:,:,ismember(psychsci_e2.times,timeIdx));%get only the times that align for all experiments
            tmpeeg_all{2}=psychsci_e2.raw_allE(:,:,ismember(psychsci_e2.times,timeIdx));%get only the times that align for all experiments
        end
        tmpbeh(2)=psychsci_e2.k;
        tmpbeh_wo_ss2(2)=psychsci_e2.k_ss4only;
        tmpbeh_only_ss2(2)=psychsci_e2.k_ss2only;
    end
    if whichExp(3)==1%if this participant did psych sci exp 3
        tmpnum=psycSci.sub_nums.e3(psycSci.uniqueIds.e3==thisId);
        load([path_psycSci,'/experiment3/',num2str(tmpnum),'_eeg_beh_psycsci_e3.mat'])
        if whichTimes==1
            tmpeeg{3}=psychsci_e3.raw_baseline(:,:,ismember(psycSci.times_baseline.e3,timeIdx));
        elseif whichTimes==2
            tmpeeg{3}=psychsci_e3.raw(:,:,ismember(psychsci_e3.times,timeIdx));%get only the times that align for all experiments
            tmpeeg_all{3}=psychsci_e3.raw_allE(:,:,ismember(psychsci_e3.times,timeIdx));%get only the times that align for all experiments
        end
        tmpbeh(3)=psychsci_e3.k;
        tmpbeh_wo_ss2(3)=psychsci_e3.k_ss4only;
        tmpbeh_only_ss2(3)=psychsci_e3.k_ss2only;
    end
    if whichExp(4)==1%if this participant did psych sci exp 4
        tmpnum=psycSci.sub_nums.e4(psycSci.uniqueIds.e4==thisId);
        load([path_psycSci,'/experiment4/',num2str(tmpnum),'_eeg_beh_psycsci_e4.mat'])
        if whichTimes==1
            tmpeeg{4}=psychsci_e4.raw_baseline(:,:,ismember(psycSci.times_baseline.e4,timeIdx));
        elseif whichTimes==2
            tmpeeg{4}=psychsci_e4.raw(:,:,ismember(psychsci_e4.times,timeIdx));%get only the times that align for all experiments
            tmpeeg_all{4}=psychsci_e4.raw_allE(:,:,ismember(psychsci_e4.times,timeIdx));%get only the times that align for all experiments
        end
        tmpbeh(4)=psychsci_e4.k;
        tmpbeh_wo_ss2(4)=psychsci_e4.k_ss4only;
        tmpbeh_only_ss2(4)=psychsci_e4.k_ss2only;
    end
    if whichExp(5)==1%if this participant did jocn exp 1
        tmpnum=jocn.sub_nums.e1(jocn.uniqueIds.e1==thisId);
        load([path_jocn,'/experiment1/',num2str(tmpnum),'_eeg_beh_jocn_e1.mat'])
        if whichTimes==1
            tmpeeg{5}=jocn_e1.raw_baseline(:,:,ismember(jocn_e2.times_baseline,timeIdx));
        elseif whichTimes==2
            tmpeeg{5}=jocn_e1.raw(:,:,ismember(jocn_e1.times,timeIdx));
            tmpeeg_all{5}=jocn_e1.raw_allE(:,:,ismember(jocn_e1.times,timeIdx));
        end
        tmpbeh(5)=jocn_e1.beh.k;
        tmpbeh_wo_ss2(5)=tmpbeh(5);%because this experiment doesn't include set size 2
    end
    if whichExp(6)==1%if this participant did jocn exp 2
        tmpnum=jocn.sub_nums.e2(jocn.uniqueIds.e2==thisId);
        load([path_jocn,'/experiment2/',num2str(tmpnum),'_eeg_beh_jocn_e2.mat'])
        if whichTimes==1
            tmpeeg{6}=jocn_e2.raw_baseline(:,:,ismember(jocn_e2.times_baseline,timeIdx));
        elseif whichTimes==2
            tmpeeg{6}=jocn_e2.raw(:,:,ismember(jocn_e2.times,timeIdx));
            tmpeeg_all{6}=jocn_e2.raw_allE(:,:,ismember(jocn_e2.times,timeIdx));
        end
        tmpbeh(6)=jocn_e2.beh.k;
        tmpbeh_wo_ss2(6)=tmpbeh(6);%because this experiment doesn't include set size 2
    end
    if whichExp(7)==1%if this participant did unpublished int exp 4
        tmpnum=unpub_int.sub_nums.e4(unpub_int.uniqueIds.e4==thisId);
        load([path_unpubint,'/experiment4/',num2str(tmpnum),'_eeg_beh_unpubint_e4.mat'])
        if whichTimes==1
            tmpeeg{7}=unpubint_e4.raw_baseline(:,:,ismember(unpubint_e4.times_baseline,timeIdx));
        elseif whichTimes==2
            tmpeeg{7}=unpubint_e4.raw(:,:,ismember(unpubint_e4.times,timeIdx));
            tmpeeg_all{7}=unpubint_e4.raw_allE(:,:,ismember(unpubint_e4.times,timeIdx));
        end
        tmpbeh(7)=unpubint_e4.beh.k;
        tmpbeh_wo_ss2(7)=tmpbeh(7);%because this experiment doesn't include set size 2
    end
    if whichExp(8)==1%if this participant did unpublished int exp 5
        tmpnum=unpub_int.sub_nums.e5(unpub_int.uniqueIds.e5==thisId);
        load([path_unpubint,'/experiment5/',num2str(tmpnum),'_eeg_beh_unpubint_e5.mat'])
        if whichTimes==1
            tmpeeg{8}=unpubint_e5.raw_baseline(:,:,ismember(unpubint_e5.times_baseline,timeIdx));
        elseif whichTimes==2
            tmpeeg{8}=unpubint_e5.raw(:,:,ismember(unpubint_e5.times,timeIdx));
            tmpeeg_all{8}=unpubint_e5.raw_allE(:,:,ismember(unpubint_e5.times,timeIdx));
        end
        tmpbeh(8)=unpubint_e5.beh.k;
        tmpbeh_wo_ss2(8)=unpubint_e5.beh.k_ss4only;%because this experiment doesn't include set size 2
        tmpbeh_only_ss2(8)=unpubint_e5.beh.k_ss2only;
    end
    if whichExp(9)==1%if this participant did relevant interruption exp 1
        tmpnum=relint.sub_nums.e1(relint.uniqueIds.e1==thisId);
        load([path_relint,'/experiment1/',num2str(tmpnum),'_eeg_beh_relint_e1.mat'])
        if whichTimes==1
            tmpeeg{9}=relint_e1.raw_baseline(:,:,ismember(relint_e1.times_baseline,timeIdx));
        elseif whichTimes==2
            tmpeeg{9}=relint_e1.raw(:,:,ismember(relint_e1.times,timeIdx));
            tmpeeg_all{9}=relint_e1.raw_allE(:,:,ismember(relint_e1.times,timeIdx));
        end
        tmpbeh(9)=relint_e1.beh.k;
        tmpbeh_wo_ss2(9)=tmpbeh(9);%because this experiment doesn't include set size 2
    end
    if whichExp(10)==1%if this participant did relevant interruption exp 2
        tmpnum=relint.sub_nums.e2(relint.uniqueIds.e2==thisId);
        load([path_relint,'/experiment2/',num2str(tmpnum),'_eeg_beh_relint_e2.mat'])
        if whichTimes==1
            tmpeeg{10}=relint_e2.raw_baseline(:,:,ismember(relint_e2.times_baseline,timeIdx));
        elseif whichTimes==2
            tmpeeg{10}=relint_e2.raw(:,:,ismember(relint_e2.times,timeIdx));
            tmpeeg_all{10}=relint_e2.raw_allE(:,:,ismember(relint_e2.times,timeIdx));
        end
        tmpbeh(10)=relint_e2.beh.k;
        tmpbeh_wo_ss2(10)=tmpbeh(10);%because this experiment doesn't include set size 2
    end
    if whichExp(11)==1%if this participant did moselle exp 1
        tmpnum=moselle.sub_nums.e1(moselle.uniqueIds.e1==thisId);
        load([path_moselle,'/experiment1/',num2str(tmpnum),'_eeg_beh_moselle_e1.mat'])
        if whichTimes==1
            tmpeeg{11}=moselle_e1.raw_baseline(:,:,ismember(moselle_e1.times_baseline,timeIdx));
        elseif whichTimes==2
            tmpeeg{11}=moselle_e1.raw(:,:,ismember(moselle_e1.times,timeIdx));
            tmpeeg_all{11}=moselle_e1.raw_allE(:,:,ismember(moselle_e1.times,timeIdx));
        end
        tmpbeh(11)=moselle_e1.k;
        tmpbeh_wo_ss2(11)=moselle_e1.k_ss4only;
        tmpbeh_only_ss2(11)=moselle_e1.k_ss2only;
    end
    if whichExp(12)==1%if this participant did moselle exp 2
        tmpnum=moselle.sub_nums.e2(moselle.uniqueIds.e2==thisId);
        load([path_moselle,'/experiment2/',num2str(tmpnum),'_eeg_beh_moselle_e2.mat'])
        if whichTimes==1
            tmpeeg{12}=moselle_e2.raw_baseline(:,:,ismember(moselle_e2.times_baseline,timeIdx));
        elseif whichTimes==2
            tmpeeg{12}=moselle_e2.raw(:,:,ismember(moselle_e2.times,timeIdx));
            tmpeeg_all{12}=moselle_e2.raw_allE(:,:,ismember(moselle_e2.times,timeIdx));
        end
        tmpbeh(12)=moselle_e2.k;
        tmpbeh_wo_ss2(12)=moselle_e2.k_ss4only;
        tmpbeh_only_ss2(12)=moselle_e2.k_ss2only;
    end
    
    tmpeeg=tmpeeg(~cellfun('isempty',tmpeeg));%remove any empty cell arrays
    tmpeeg_all=tmpeeg_all(~cellfun('isempty',tmpeeg_all));%remove any empty cell arrays
    %%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %for everyone, I also want to average over all studies to get one fc
    %matrix per person -- This is what happens below
    %%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %average the behavioral results and save to fc_mat.beh
    fc_mat.beh(s)=nanmean(tmpbeh);
    fc_mat.beh_wo_ss2(s)=nanmean(tmpbeh_wo_ss2); 
    fc_mat.beh_only_ss2(s)=nanmean(tmpbeh_only_ss2);
    
    %get one eeg matrix that is electrodes x times (by averaging over the trial info from all experiments)
    %only overlapping electrodes 
    tmpeeg1=cat(1,tmpeeg{:});
    fc_mat.nTrials(s)=size(tmpeeg1,1);%save the number of trials per participant
    tmpeeg=cat(1,tmpeeg{:});
    tmpeeg=squeeze(nanmean(tmpeeg,1));%concatenate all of the trials and then average over them
    
    %save the raw amplitude of the relevant trials
    fc_mat.amplitude.raw(s,:,:)=tmpeeg;
    fc_mat.amplitude.raw_aveTime(s,:)=nanmean(tmpeeg,2);
    
    %all electrodes
    tmpeeg1_all=cat(1,tmpeeg_all{:});
    tmpeeg_all=cat(1,tmpeeg_all{:});
    tmpeeg_all=squeeze(nanmean(tmpeeg_all,1));%concatenate all of the trials and then average over them
    
    %save contra - ipsi
    %all electrodes
    tmp_cmi=nan(1,size(contra_ipsi_elecs,1));
    tmp_cmi_laterTimes=nan(1,size(contra_ipsi_elecs,1));
    for pair=1:size(contra_ipsi_elecs,1)%loop through the contra ipsi pairs and subtract the difference
        tmp_cmi(pair)=nanmean(nanmean(tmpeeg1(:,ismember(raw_order,contra_ipsi_elecs{pair,1}),:)-tmpeeg1(:,ismember(raw_order,contra_ipsi_elecs{pair,2}),:),1),3);
        tmp_cmi_laterTimes(pair)=nanmean(nanmean(tmpeeg1(:,ismember(raw_order,contra_ipsi_elecs{pair,1}),ismember(timeIdx,400:1000))-tmpeeg1(:,ismember(raw_order,contra_ipsi_elecs{pair,2}),ismember(timeIdx,400:1000)),1),3);
    end
    fc_mat.contraMipsi(s,:)=tmp_cmi;
    fc_mat.contraMipsi_400_to_1000(s,:)=tmp_cmi;
    
    %now, make the functional connectivity matrix from the eeg data 
    %for only overlapping electrodes
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                fc_mat.eeg(s,ed1,ed2)=18;
            else
                fc_mat.eeg(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,:)',tmpeeg(ed2,:)','type','Pearson'));
            end
        end
    end
    
    if whichTimes==2 %if this the retention interval window, then make an fc_mat for all electrodes
        %now, make the functional connectivity matrix from the eeg data
        %for all electrodes
        for ed1=1:30
            for ed2=1:30
                if ed1==ed2
                    fc_mat.eeg_allE(s,ed1,ed2)=18;
                else
                    fc_mat.eeg_allE(s,ed1,ed2)=atanh(corr(tmpeeg_all(ed1,:)',tmpeeg_all(ed2,:)','type','Pearson'));
                end
            end
        end
    end
    fc_mat.experiment(s,1:sum(whichExp))=find(whichExp);%save which experiment they participated in
    
    toc
end

%save the functional connectivity data 
save([savePath,saveName,'.mat'],'fc_mat')




