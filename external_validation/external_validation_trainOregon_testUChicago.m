%% This script externally validates the Oregon-site model
clear;clc

%% Some settings
dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/data/';%location where the data is saved
nIter=10000;%used to calculate shuffled r and mse
offset=1;%if you want to train and test the model from stimulus offset (0=no [0-1000 ms]; 1=yes [offset-offset+800ms]

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
train_mats  = oregon_fc;                     % training data (n_node x n_node x n_sub symmetrical connectivity matrices)
behav       = oregon_behav;                  % n_sub x 1 vector of behavior
n_sub       = size(oregon_fc,3);             % total number of training subjects
n_train_sub = n_sub-1;                       % number of training subjects in each round of leave-one-out

% Validation data
validation_mats  = uchicago_fc;              % validation data (n_node x n_node x n_validation_sub symmetrical connectivity matrices)
validation_behav = uchicago_beh;             % n_validation_sub x 1 vector of behavior
n_validation_sub = size(validation_mats,3);  % total number of validation subjects

aa     = ones(n_node,n_node);
aa_upp = triu(aa,1);
upp_id = find(aa_upp);   % indices of edges in the upper triangular of an n_node x n_node matrix
n_edge = length(upp_id); % total number of edges

%% Feature selection: use top 10% of edges
load([dir,'oregon_chicago_top10percedges']);
importantEdges=triu(sigEdges_overlap.color_shape,1);%make sure that only the upper half of the importantedges are included (so that we don't double count a specific edge

pos_mask = importantEdges==1;
neg_mask = importantEdges==-1;

pos_overlap = pos_mask;
neg_overlap = neg_mask;
%% Model building
% sum edges for all subjects in the training set
train_pos_sum = zeros(n_sub,1);
train_neg_sum = zeros(n_sub,1);
    
for k = 1:n_sub
    train_pos_sum(k) = nansum(nansum(pos_overlap.*train_mats(:,:,k)));
    train_neg_sum(k) = nansum(nansum(neg_overlap.*train_mats(:,:,k)));
end
    
% build model with training data
robGLM_fit = robustfit(train_pos_sum-train_neg_sum,behav);
    
%% Generate predictions for validation set
pred_glm = zeros(n_validation_sub,1);

validation_pos_sum = zeros(n_validation_sub,1);
validation_neg_sum = zeros(n_validation_sub,1);

for vs = 1:n_validation_sub
    validation_pos_sum(vs) = sum(sum(pos_overlap.*validation_mats(:,:,vs)));
    validation_neg_sum(vs) = sum(sum(neg_overlap.*validation_mats(:,:,vs)));

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

%% Plot results
figure 
scatter(validation_behav, pred_glm,200,'filled','MarkerFaceAlpha',.4)
xlabel('Observed')
ylabel('Predicted')
set(gca,'FontSize',30,'LineWidth',5)
title({'Predicting Chicago K scores from Oregon model',['glm r = ' num2str(round(r_glm*100)/100) ', p = ' num2str(p_glm),],[' mse = ',num2str(mse),', p < ',num2str(mse_p_para)]})

%% Plot distribution of shuffled MSE
lim_bar=nIter/20;

figure; hold on 
bar(mse,lim_bar,'FaceColor',[0 0 0],'BarWidth',0.003,'FaceAlpha',1,'EdgeAlpha',0);
histogram(shuff_mse,'FaceColor',[.65 .65 .65],'EdgeAlpha',0)%histogram of shuffled results
set(gca,'FontSize',20,'LineWidth',2)
box off
ylabel('Count')
xlabel('predicted r squared')
legend('predicted r squared','null distribution','Location','northwest');legend('boxoff');
title({'Train Oregon, test Chicago','Mean square error',['mse=',num2str(mse),', p<',num2str(mse_p_para)]})
%% unique colors 
uniqColors=distinguishable_colors(12);
%% Plot GLM dots separated by experiment
pred_glm = zeros(n_validation_sub,1);
x=1:5;

validation_pos_sum = zeros(n_validation_sub,1);
validation_neg_sum = zeros(n_validation_sub,1);

for vs = 1:n_validation_sub
    validation_pos_sum(vs) = sum(sum(pos_overlap.*validation_mats(:,:,vs)));
    validation_neg_sum(vs) = sum(sum(neg_overlap.*validation_mats(:,:,vs)));

    pred_glm(vs) = robGLM_fit(1) + robGLM_fit(2)*(validation_pos_sum(vs)-validation_neg_sum(vs));
end

%make plot
figure; set(gcf,'Position',[57 1000 1000 800]); hold on

scatter(validation_behav, pred_glm,200,[0 0 0],'filled','MarkerFaceAlpha',.5,'LineWidth',5)%uniqColors(e,:),'filled','MarkerFaceAlpha',.6)

% correlate predicted and observed behavior
[r_glm, p_glm] = corr(validation_behav, pred_glm,'type','spearman');
mse_across_studies=nan(12,1);
for e=1:12
    exp_ind=logical(sum(ismember(uchicago.fc_mat.experiment,e),2));
    mse_across_studies(e)=mean((validation_behav(exp_ind)-pred_glm(exp_ind)).^2);%calculate mean square error for each study
    
    [y_pred,delta_pred]=polyfit(validation_behav(exp_ind), pred_glm(exp_ind),1);
    y_fit=y_pred(1)*x+y_pred(2);
    
    p1=plot(x,y_fit,'LineWidth',5,'Color',uniqColors(e,:));p1.Color(4)=0.35;
end

%regression line for all points
[y_pred,delta_pred]=polyfit(validation_behav, pred_glm,1);
y_fit_all=y_pred(1)*x+y_pred(2);
plot(x,y_fit_all,'LineWidth',10,'Color',[0 0 0]);

if offset==0
    title({'Train Oregon, test UChicago',['glm r = ' num2str(round(r_glm*100)/100) ', p = ' num2str(p_glm),],[' mse = ',num2str(mse),', p = ',num2str(mse_p_para)]})
elseif offset==1
    title({'Train Oregon, test UChicago','offset-based model',['glm r = ' num2str(round(r_glm*100)/100) ', p = ' num2str(p_glm),],[' mse = ',num2str(mse),', p = ',num2str(mse_p_para)]})
end

xlabel('Observed')
ylabel('Predicted')
set(gca,'FontSize',45,'LineWidth',2)
ylim([1 5]);xlim([1 5])


