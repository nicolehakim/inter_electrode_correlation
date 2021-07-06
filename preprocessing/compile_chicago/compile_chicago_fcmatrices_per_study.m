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

importantTimes=0:1000;
%% Path settings
path_psycSci='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_PsychSci/data/';
path_jocn='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_JoCN/data/';
path_unpubint='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_unpublished_int/data/';
path_relint='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Hakim_RelInt/data/';
path_moselle='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/Toby_Moselle/data/';
%% Load the subject numbers and unique ids for each experiment 
%psych science
load([path_psycSci,'PsycSci_uniqueIds.mat']);
uniqueSubs.psycSci=unique([psycSci.sub_nums.e1,psycSci.sub_nums.e2,psycSci.sub_nums.e3,psycSci.sub_nums.e4]);%unique(psycSci.subs.uniqueIds(~isnan(psycSci.subs.uniqueIds)));
% psycSci.sub_nums.e4(27)=37;%this subject number was incorrect in the original experiment file. this will be fixed once everything re-runs

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
uniqueSubs.all=unique([uniqueSubs.psycSci,uniqueSubs.jocn,uniqueSubs.unpub,uniqueSubs.relint,uniqueSubs.moselle]);
nUniqueSubs=length(uniqueSubs.all);

%place to save this data
savePath='/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/';

%% Order of chicago electrodes
raw_order={'F3','F4','C3','C4','P7','P8','P3','P4',...
    'PO7','PO8','PO3','PO4','O1','O2','Fz','Cz','Pz'};
% save the raw contra - ipsi electrodes
contra_ipsi_elecs={'PO3','PO4';'F3','F4';'C3','C4';'P3','P4';'O1','O2';'PO7','PO8';'P7','P8'};
nElec=17; %number of overlapping electrodes between uchicago and oregon

%% Save fc_mat for the psych sci exp 1 dataset
psychSci1.expNum=1; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)

tmpSubs=psycSci.sub_nums.e1;tmpSubs=tmpSubs(~isnan(tmpSubs));
nSubs=length(tmpSubs);

psychSci1.uniqID=psycSci.uniqueIds.e1;
psychSci1.elecOrder=raw_order; 
psychSci1.contraIpsiOrder=contra_ipsi_elecs;

psychSci1.fc_mat=nan(nSubs,nElec,nElec);
psychSci1.beh=nan(nSubs,1);
psychSci1.whichExps=nan(nSubs,12);

fprintf('Psych Sci 1 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_psycSci,'experiment1/',num2str(tmpSubs(s)),'_eeg_beh_psycsci_e1.mat'])
    thisId=psychSci1.uniqID(s);
    psychSci1.beh(s)=psychsci_e1.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(psychsci_e1.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                psychSci1.fc_mat(s,ed1,ed2)=18; %because it will be Inf otherwise
            else
                psychSci1.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(psychsci_e1.times,importantTimes))',tmpeeg(ed2,ismember(psychsci_e1.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    psychSci1.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp1_psychSci1.mat'],'psychSci1')

%% Save fc_mat for the psych sci exp 2 dataset
psychSci2.expNum=2; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)

tmpSubs=psycSci.sub_nums.e2;tmpSubs=tmpSubs(~isnan(tmpSubs));
nSubs=length(tmpSubs);

psychSci2.uniqID=psycSci.uniqueIds.e2;
psychSci2.elecOrder=raw_order; 
psychSci2.contraIpsiOrder=contra_ipsi_elecs;

psychSci2.fc_mat=nan(nSubs,nElec,nElec);
psychSci2.beh=nan(nSubs,1);
psychSci2.whichExps=nan(nSubs,12);

fprintf('Psych Sci 2 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_psycSci,'experiment2/',num2str(tmpSubs(s)),'_eeg_beh_psycsci_e2.mat'])
    thisId=psychSci2.uniqID(s);
    psychSci2.beh(s)=psychsci_e2.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(psychsci_e2.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                psychSci2.fc_mat(s,ed1,ed2)=18; 
            else
                psychSci2.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(psychsci_e2.times,importantTimes))',tmpeeg(ed2,ismember(psychsci_e2.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    psychSci2.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp2_psychSci2.mat'],'psychSci2')

%% Save fc_mat for the psych sci exp 3 dataset
psychSci3.expNum=3; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)

tmpSubs=psycSci.sub_nums.e3;tmpSubs=tmpSubs(~isnan(tmpSubs));
nSubs=length(tmpSubs);

psychSci3.uniqID=psycSci.uniqueIds.e3;
psychSci3.elecOrder=raw_order; 
psychSci3.contraIpsiOrder=contra_ipsi_elecs;

psychSci3.fc_mat=nan(nSubs,nElec,nElec);
psychSci3.beh=nan(nSubs,1);
psychSci3.whichExps=nan(nSubs,12);

fprintf('Psych Sci 3 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_psycSci,'experiment3/',num2str(tmpSubs(s)),'_eeg_beh_psycsci_e3.mat'])
    thisId=psychSci3.uniqID(s);
    psychSci3.beh(s)=psychsci_e3.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(psychsci_e3.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                psychSci3.fc_mat(s,ed1,ed2)=18; 
            else
                psychSci3.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(psychsci_e3.times,importantTimes))',tmpeeg(ed2,ismember(psychsci_e3.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    psychSci3.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp3_psychSci3.mat'],'psychSci3')

%% Save fc_mat for the psych sci exp 4 dataset
psychSci4.expNum=4; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)

tmpSubs=psycSci.sub_nums.e4;tmpSubs=tmpSubs(~isnan(tmpSubs));
nSubs=length(tmpSubs);

psychSci4.uniqID=psycSci.uniqueIds.e4;
psychSci4.elecOrder=raw_order; 
psychSci4.contraIpsiOrder=contra_ipsi_elecs;

psychSci4.fc_mat=nan(nSubs,nElec,nElec);
psychSci4.beh=nan(nSubs,1);
psychSci4.whichExps=nan(nSubs,12);

fprintf('Psych Sci 4 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_psycSci,'experiment4/',num2str(tmpSubs(s)),'_eeg_beh_psycsci_e4.mat'])
    thisId=psychSci4.uniqID(s);
    psychSci4.beh(s)=psychsci_e4.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(psychsci_e4.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                psychSci4.fc_mat(s,ed1,ed2)=18; 
            else
                psychSci4.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(psychsci_e4.times,importantTimes))',tmpeeg(ed2,ismember(psychsci_e4.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    psychSci4.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp4_psychSci4.mat'],'psychSci4')

%% Save fc_mat for the jocn exp 1
jocn1.expNum=5; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)

tmpSubs=jocn.sub_nums.e1;
nSubs=length(tmpSubs);
jocn1.uniqID=jocn.uniqueIds.e1;

jocn1.elecOrder=raw_order; 
jocn1.contraIpsiOrder=contra_ipsi_elecs;

jocn1.fc_mat=nan(nSubs,nElec,nElec);
jocn1.beh=nan(nSubs,1);
jocn1.whichExps=nan(nSubs,12);

fprintf('JoCN 1 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_jocn,'experiment1/',num2str(tmpSubs(s)),'_eeg_beh_jocn_e1.mat'])
    thisId=jocn1.uniqID(s);
    jocn1.beh(s)=jocn_e1.beh.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(jocn_e1.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                jocn1.fc_mat(s,ed1,ed2)=18; 
            else
                jocn1.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(jocn_e1.times,importantTimes))',tmpeeg(ed2,ismember(jocn_e1.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    jocn1.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp5_jocn1.mat'],'jocn1')

%% Save fc_mat for the jocn exp 2
jocn2.expNum=6; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)
tmpSubs=jocn.sub_nums.e2;
nSubs=length(tmpSubs);
jocn2.uniqID=jocn.uniqueIds.e2;

jocn2.elecOrder=raw_order; 
jocn2.contraIpsiOrder=contra_ipsi_elecs;
jocn2.fc_mat=nan(nSubs,nElec,nElec);
jocn2.beh=nan(nSubs,1);
jocn2.whichExps=nan(nSubs,12);

fprintf('JoCN 2 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_jocn,'experiment2/',num2str(tmpSubs(s)),'_eeg_beh_jocn_e2.mat'])
    thisId=jocn2.uniqID(s);
    jocn2.beh(s)=jocn_e2.beh.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(jocn_e2.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                jocn2.fc_mat(s,ed1,ed2)=18; 
            else
                jocn2.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(jocn_e2.times,importantTimes))',tmpeeg(ed2,ismember(jocn_e2.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    jocn2.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp6_jocn2.mat'],'jocn2')

%% Save fc_mat for the unpublished int 4
unpub1.expNum=7; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)
tmpSubs=unpub_int.sub_nums.e4;
nSubs=length(tmpSubs);
unpub1.uniqID=unpub_int.uniqueIds.e4;

unpub1.elecOrder=raw_order; 
unpub1.contraIpsiOrder=contra_ipsi_elecs;
unpub1.fc_mat=nan(nSubs,nElec,nElec);
unpub1.beh=nan(nSubs,1);
unpub1.whichExps=nan(nSubs,12);

fprintf('Unpublished int 1 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_unpubint,'experiment4/',num2str(tmpSubs(s)),'_eeg_beh_unpubint_e4.mat'])
    thisId=unpub1.uniqID(s);
    unpub1.beh(s)=unpubint_e4.beh.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(unpubint_e4.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                unpub1.fc_mat(s,ed1,ed2)=18; 
            else
                unpub1.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(unpubint_e4.times,importantTimes))',tmpeeg(ed2,ismember(unpubint_e4.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    unpub1.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp7_unpub1.mat'],'unpub1')

%% Save fc_mat for the unpublished int 5
unpub2.expNum=8; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)
tmpSubs=unpub_int.sub_nums.e5;
nSubs=length(tmpSubs);
unpub2.uniqID=unpub_int.uniqueIds.e5;

unpub2.elecOrder=raw_order; 
unpub2.contraIpsiOrder=contra_ipsi_elecs;
unpub2.fc_mat=nan(nSubs,nElec,nElec);
unpub2.beh=nan(nSubs,1);
unpub2.whichExps=nan(nSubs,12);

fprintf('Unpublished int 2 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_unpubint,'experiment5/',num2str(tmpSubs(s)),'_eeg_beh_unpubint_e5.mat'])
    unpub2.beh(s)=unpubint_e5.beh.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(unpubint_e5.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                unpub2.fc_mat(s,ed1,ed2)=18; 
            else
                unpub2.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(unpubint_e5.times,importantTimes))',tmpeeg(ed2,ismember(unpubint_e5.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    thisId=unpub2.uniqID(s);
    unpub2.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp8_unpub2.mat'],'unpub2')

%% Save fc_mat for the relevant int exp 1
relint1.expNum=9; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)
tmpSubs=relint.sub_nums.e1;
nSubs=length(tmpSubs);
relint1.uniqID=relint.uniqueIds.e1;

relint1.elecOrder=raw_order; 
relint1.contraIpsiOrder=contra_ipsi_elecs;
relint1.fc_mat=nan(nSubs,nElec,nElec);
relint1.beh=nan(nSubs,1);
relint1.whichExps=nan(nSubs,12);

fprintf('Relevant interruption 1 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_relint,'experiment1/',num2str(tmpSubs(s)),'_eeg_beh_relint_e1.mat'])
    relint1.beh(s)=relint_e1.beh.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(relint_e1.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                relint1.fc_mat(s,ed1,ed2)=18; 
            else
                relint1.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(relint_e1.times,importantTimes))',tmpeeg(ed2,ismember(relint_e1.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    thisId=relint1.uniqID(s);
    relint1.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp9_relint1.mat'],'relint1')

%% Save fc_mat for the relevant int exp 2
relint2.expNum=10; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)
tmpSubs=relint.sub_nums.e2;
nSubs=length(tmpSubs);
relint2.uniqID=relint.uniqueIds.e2;

relint2.elecOrder=raw_order; 
relint2.contraIpsiOrder=contra_ipsi_elecs;
relint2.fc_mat=nan(nSubs,nElec,nElec);
relint2.beh=nan(nSubs,1);
relint2.whichExps=nan(nSubs,12);

fprintf('Relevant interruption 2 \n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_relint,'experiment2/',num2str(tmpSubs(s)),'_eeg_beh_relint_e2.mat'])
    relint2.beh(s)=relint_e2.beh.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(relint_e2.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                relint2.fc_mat(s,ed1,ed2)=18; 
            else
                relint2.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(relint_e2.times,importantTimes))',tmpeeg(ed2,ismember(relint_e2.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    thisId=relint2.uniqID(s);
    relint2.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp10_relint2.mat'],'relint2')

%% Save fc_mat for the Moselle exp 1
moselle1.expNum=11; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)
tmpSubs=moselle.sub_nums.e1;
nSubs=length(tmpSubs);
moselle1.uniqID=moselle.uniqueIds.e1;

moselle1.elecOrder=raw_order; 
moselle1.contraIpsiOrder=contra_ipsi_elecs;
moselle1.fc_mat=nan(nSubs,nElec,nElec);
moselle1.beh=nan(nSubs,1);
moselle1.whichExps=nan(nSubs,12);

fprintf('Moselle 1\n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_moselle,'experiment1/',num2str(tmpSubs(s)),'_eeg_beh_moselle_e1.mat'])
    moselle1.beh(s)=moselle_e1.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(moselle_e1.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                moselle1.fc_mat(s,ed1,ed2)=18; 
            else
                moselle1.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(moselle_e1.times,importantTimes))',tmpeeg(ed2,ismember(moselle_e1.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    thisId=moselle1.uniqID(s);
    moselle1.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp11_moselle1.mat'],'moselle1')

%% Save fc_mat for the Moselle exp 2
moselle2.expNum=12; %the number of this experiment (given arbitarily, so I can keep track of each experiment using a number, rather than its name)
tmpSubs=moselle.sub_nums.e2;
nSubs=length(tmpSubs);
moselle2.uniqID=moselle.uniqueIds.e2;

moselle2.elecOrder=raw_order; 
moselle2.contraIpsiOrder=contra_ipsi_elecs;
moselle2.fc_mat=nan(nSubs,nElec,nElec);
moselle2.beh=nan(nSubs,1);
moselle2.whichExps=nan(nSubs,12);

fprintf('Moselle 2\n')
for s=1:nSubs
    fprintf(['Subject ',num2str(s),'\n'])
    load([path_moselle,'experiment2/',num2str(tmpSubs(s)),'_eeg_beh_moselle_e2.mat'])
    moselle2.beh(s)=moselle_e2.k;
    
    %average over trials 
    tmpeeg=squeeze(mean(moselle_e2.raw,1));
    %calculate connectivity matrix 
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1==ed2
                moselle2.fc_mat(s,ed1,ed2)=18; 
            else
                moselle2.fc_mat(s,ed1,ed2)=atanh(corr(tmpeeg(ed1,ismember(moselle_e2.times,importantTimes))',tmpeeg(ed2,ismember(moselle_e2.times,importantTimes))','type','Pearson'));
            end
        end
    end
    
    %save the experiment that this participant did 
    thisId=moselle2.uniqID(s);
    moselle2.whichExps(s,:)=[ismember(thisId,psycSci.uniqueIds.e1),...  %Psych sci exp 1
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
end

save([savePath,'fc_mat_exp12_moselle2.mat'],'moselle2')