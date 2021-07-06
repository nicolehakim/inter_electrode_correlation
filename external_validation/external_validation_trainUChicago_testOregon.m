%% This script externally validates the Oregon-site model
clear;clc

%% Some settings
dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/data/';%location where the data is saved
nIter=10000;%used to calculate shuffled r and mse
offset=0;%if you want to train and test the model from stimulus offset (0=no [0-1000 ms]; 1=yes [offset-offset+800ms]

%% Load & correctly configure the data

if offset==0%include timepoints 0-1000 ms for all studies 
    %load the Oregon-site data
    load([dir,'oregon_site_fc_mat.mat']);
    %load the Chicago-site data
    uchicago=load([dir,'chicago_site_fc_mat.mat']);
elseif offset==1%yes, only include data from stimulus offset
    %load the Oregon-site data
    bs=load([dir,'oregon_site_fc_mat_from_offset.mat']);
    %load the Chicago-site data
    uchicago=load([dir,'chicago_site_fc_mat_from_offset.mat']); 
end

%correctly configure the Oregon-site data
oregon_fc    = mean(cat(4,bs.fc_mat.overlappingElecs.RI.color,bs.fc_mat.overlappingElecs.RI.shape),4);
oregon_fc    = permute(oregon_fc,[2,3,1]);
oregon_behav = bs.fc_mat.k.ave;

%correctly configure the Chicago-site data
uchicago_fc  = uchicago.fc_mat.eeg;
uchicago_fc    = permute(uchicago_fc,[2,3,1]);
uchicago_beh = uchicago.fc_mat.beh; % Make sure subject and fc data indices align

%% Define the training and testing datasets
n_node  = size(uchicago_fc,1);   % number of nodes

% Training data 
train_mats  = uchicago_fc;                      % training data (n_node x n_node x n_sub symmetrical connectivity matrices)
behav       = uchicago_beh;                     % n_sub x 1 vector of behavior
n_sub       = size(uchicago_fc,3);              % total number of training subjects
n_train_sub = n_sub-1;                          % number of training subjects in each round of leave-one-out

% Validation data
validation_mats  = oregon_fc;                   % validation data (n_node x n_node x n_validation_sub symmetrical connectivity matrices)
validation_behav = oregon_behav;                % n_validation_sub x 1 vector of behavior
n_validation_sub = size(validation_mats,3);     % total number of validation subjects

aa     = ones(n_node,n_node);
aa_upp = triu(aa,1);
upp_id = find(aa_upp);   % indices of edges in the upper triangular of an n_node x n_node matrix
n_edge = length(upp_id); % total number of edges

%% Feature selection: use top 10% of edges
load([dir,'oregon_chicago_top10percedges']);
importantEdges=triu(sigEdges.chicago,1); %make sure that only the upper half of the importantedges are included (so that we don't double count a specific edge

pos_mask = importantEdges==1;
neg_mask = importantEdges==-1;

pos_overlap = pos_mask;
neg_overlap = neg_mask;
%% Model building
% sum edges for all subjects in the training set
train_pos_sum = zeros(n_sub,1);
train_neg_sum = zeros(n_sub,1);
    
for k = 1:n_sub
    train_pos_sum(k) = sum(sum(pos_overlap.*train_mats(:,:,k)));
    train_neg_sum(k) = sum(sum(neg_overlap.*train_mats(:,:,k)));
end
    
% build model with training data
robGLM_fit = robustfit(train_pos_sum-train_neg_sum,behav);
    
%% Generate predictions for validation set
pred_glm = zeros(n_validation_sub,1);

validation_pos_sum = zeros(n_validation_sub,1);
validation_neg_sum = zeros(n_validation_sub,1);

for vs = 1:n_validation_sub
    validation_pos_sum(vs) = nansum(nansum(pos_overlap.*validation_mats(:,:,vs)));
    validation_neg_sum(vs) = nansum(nansum(neg_overlap.*validation_mats(:,:,vs)));

    pred_glm(vs) = robGLM_fit(1) + robGLM_fit(2)*(validation_pos_sum(vs)-validation_neg_sum(vs));
end

% correlate predicted and observed behavior
[r_glm, p_glm] = corr(validation_behav, pred_glm,'Type','Spearman');

%% create null distribution & calculate non-parametric p-value for r
z_shuff=nan(1,nIter);
for ii=1:nIter
    z_shuff(ii)=atanh(corr(shuffle(validation_behav), pred_glm,'Type','Spearman'));
end
r_shuff=tanh(z_shuff);
p_para=(1+sum(r_shuff>=r_glm))/(1+nIter);

%% calculate MSE & shuffled MSE & calculate non-parametric p-value for mse
mse=immse(validation_behav, pred_glm);

shuff_mse=nan(1,nIter);
for ii=1:nIter
    shuff_mse(ii) = immse(validation_behav,shuffle(pred_glm));
end

mse_p_para=(1+sum(shuff_mse(:)<=mse))/(1+nIter);

%% plot results
figure; set(gcf,'Position',[57 1000 1000 800]); hold on

scatter(validation_behav, pred_glm,200,[0 0 0],'filled','MarkerFaceAlpha',.5)
x=1:5;
[y_pred,delta_pred]=polyfit(validation_behav, pred_glm,1);
y_fit_all=y_pred(1)*x+y_pred(2);
plot(x,y_fit_all,'LineWidth',5,'Color',[0 0 0]);
if offset==0
    title({'Predicting Oregon K scores from Chicago model',['glm r = ' num2str(round(r_glm*100)/100) ', p = ' num2str(p_glm),],[' mse = ',num2str(mse),', p = ',num2str(mse_p_para)]})
elseif offset==1
    title({'Predicting Oregon K scores from Chicago model','offset-based model',['glm r = ' num2str(round(r_glm*100)/100) ', p = ' num2str(p_glm),],[' mse = ',num2str(mse),', p = ',num2str(mse_p_para)]})
end
xlabel('Observed')
ylabel('Predicted')
set(gca,'FontSize',45,'LineWidth',2)
ylim([1 5]);xlim([1 5])


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Behavioral cross-correlation (supplemental analysis)
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/BigStudy/Data/all_tasks_dat.mat')
temp=load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/BigStudy/Data/BigStudy_filtered_erp_singleTrial.mat');%need to load this for the EEG subject numbers
temp2=load('/Volumes/nicole passport/BigStudy/K_cda_changedetection.mat');
eeg_subs=temp2.cd_cda.subNum(ismember(temp2.cd_cda.subNum,temp.data.subs));
beh_sub_idx=ismember(dat.subjects,eeg_subs);
eeg_sub_idx=ismember(eeg_subs,dat.subjects);
pred_dat=pred_glm(eeg_sub_idx);
actual_dat=validation_behav(eeg_sub_idx);

% correlation between gF and predicted K 
[r,p]=corr(validation_behav(eeg_sub_idx),dat.gF(beh_sub_idx))

task_names={'Antisaccade','ospan','rspan','symmspan','span_ave','FreeRecall',...
    'FreeRecall_repeats','PairedAssociates','PictureSource','Ravens','Cattell',...
    'NumSeries','gF','sd_col1','sd_col6','sd_mot1','sd_mot6','sd_orient1','sd_orient6',...
    'sd_shape1','sd_shape6','sd_space1','sd_space6','sd_1','sd_6'};

r_pred_other=nan(length(task_names),1);
p_pred_other=nan(length(task_names),1);

r_actual_other=nan(length(task_names),1);
p_actual_other=nan(length(task_names),1);

for ii=1:length(task_names)
    beh_var=dat.(task_names{ii});
    
    [r_pred_other(ii),p_pred_other(ii)]=corr(beh_var(beh_sub_idx),pred_dat,'type','Pearson');
    [r_actual_other(ii),p_actual_other(ii)]=corr(beh_var(beh_sub_idx),eeg_dat,'type','Pearson');
end

%% Plot correlation between predicted K and gF
figure; set(gcf,'Position',[57 1000 1000 800]); hold on

scatter(pred_dat,dat.gF(beh_sub_idx),200,[0 0 0],'filled','MarkerFaceAlpha',.5)
x=2:3.5;
[y_pred,delta_pred]=polyfit(pred_dat,dat.gF(beh_sub_idx),1);
y_fit_all=y_pred(1)*x+y_pred(2);
plot(x,y_fit_all,'LineWidth',5,'Color',[0 0 0]);

ylabel('gF')
xlabel('Predicted K')
set(gca,'FontSize',45,'LineWidth',2)
[r,p]=corr(pred_dat,dat.gF(beh_sub_idx));
