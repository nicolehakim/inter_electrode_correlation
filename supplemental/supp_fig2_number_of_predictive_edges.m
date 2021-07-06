%% This script varies the number of predictive edges included and runs the model 
clear all

%% Settings
percentage_range=.05:.05:1;%the percentage of electrodes that we will include in this analysis

dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg';

%load the oregon data
load([dir,'/data/compiled/oregon_site_fc_mat.mat']);
oregon_fc    = mean(cat(4,bs.fc_mat.overlappingElecs.RI.color,bs.fc_mat.overlappingElecs.RI.shape),4);
oregon_fc    = permute(oregon_fc,[2,3,1]);
oregon_behav = bs.fc_mat.k.ave;

%load the chicago data 
uchicago=load([dir,'/data/compiled/chicago_site_fc_mat.mat']);
uchicago_fc  = uchicago.fc_mat.eeg;
uchicago_fc    = permute(uchicago_fc,[2,3,1]);
uchicago_beh = uchicago.fc_mat.beh; % Make sure subject and fc data indices align

n_node  = size(uchicago_fc,1);%17;   % number of nodes

%% Loop through the percentages
%pre-allocate matrices to save the data 
r_pp=nan(length(percentage_range),2);% 1=train oregon, test chicago; 2=train chicago, test oregon
p_pp=nan(length(percentage_range),2);
mse_pp=nan(length(percentage_range),2);

for pp=1:length(percentage_range)
    fprintf(['Percent number ',num2str(pp),' out of ',num2str(length(percentage_range)),'\n'])
    
    %load the edges
    load(['/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/oregon_chicago_top',num2str(percentage_range(pp)*100),'percedges.mat']);
    
    for exp=1:2%train Oregon test Chicago & vice versa
        
        if exp==1 %train oregon, test chicago
            importantEdges=sigEdges_overlap.color_shape;
            
            train_mats=oregon_fc;
            behav=oregon_behav;
            n_sub=length(oregon_behav);
            
            validation_mats=uchicago_fc; 
            validation_behav=uchicago_beh;  
            n_validation_sub=length(uchicago_beh);
            
        else      %train chicago, test oregon
            importantEdges=sigEdges.chicago;
            
            train_mats=uchicago_fc;
            behav=uchicago_beh;
            n_sub=length(uchicago_beh); 
            
            validation_mats=oregon_fc;
            validation_behav = oregon_behav;  
            n_validation_sub=length(oregon_behav);
        end
        
        pos_overlap = importantEdges==1;
        neg_overlap = importantEdges==-1;
            
        %build the model
        train_pos_sum = zeros(n_sub,1);
        train_neg_sum = zeros(n_sub,1);
        
        for k = 1:n_sub
            train_pos_sum(k) = sum(sum(pos_overlap.*train_mats(:,:,k)));
            train_neg_sum(k) = sum(sum(neg_overlap.*train_mats(:,:,k)));
        end
        robGLM_fit = robustfit(train_pos_sum-train_neg_sum,behav);
        
        % Generate predictions for validation set
        pred_glm = zeros(n_validation_sub,1);
        
        validation_pos_sum = zeros(n_validation_sub,1);
        validation_neg_sum = zeros(n_validation_sub,1);
        
        for vs = 1:n_validation_sub
            validation_pos_sum(vs) = sum(sum(pos_overlap.*validation_mats(:,:,vs)));
            validation_neg_sum(vs) = sum(sum(neg_overlap.*validation_mats(:,:,vs)));
            
            pred_glm(vs) = robGLM_fit(1) + robGLM_fit(2)*(validation_pos_sum(vs)-validation_neg_sum(vs));
        end
        
        % correlate predicted and observed behavior
        [r_pp(pp,exp), p_pp(pp,exp)] = corr(validation_behav, pred_glm,'Type','Spearman');
        mse_pp(pp,exp)=immse(validation_behav, pred_glm);
        
    end
end

%% Plot results 

figure; set(gcf,'Position',[57 1000 1000 1200]);

% R value
subplot(3,1,1); hold on
plot(percentage_range*100,r_pp(:,1),'Color',[34,139,34]/256,'LineWidth',5)
plot(percentage_range*100,r_pp(:,2),'Color',[128,0,0]/256,'LineWidth',5)
ylim([0 0.4])
set(gca,'XTick',percentage_range*100,'XTickLabel',{'5','','15','','25','','35','','45','','55','','65','','75','','85','','95',''})
set(gca,'TickDir','out','FontSize',25,'LineWidth',2);
title('r value')
legend('Train Oregon test Chicago','Train Chicago test Oregon')
legend boxoff
box off 
ylabel('r value')

% p value
subplot(3,1,2); hold on
plot(percentage_range*100,p_pp(:,1),'Color',[34,139,34]/256,'LineWidth',5)
plot(percentage_range*100,p_pp(:,2),'Color',[128,0,0]/256,'LineWidth',5)
ylim([0 0.02])
set(gca,'XTick',percentage_range*100,'XTickLabel',{'5','','15','','25','','35','','45','','55','','65','','75','','85','','95',''})
title('p value')
set(gca,'TickDir','out','FontSize',25,'LineWidth',2)
box off 
ylabel('p value')

% Mean square error 
subplot(3,1,3); hold on
plot(percentage_range*100,mse_pp(:,1),'Color',[34,139,34]/256,'LineWidth',5)
plot(percentage_range*100,mse_pp(:,2),'Color',[128,0,0]/256,'LineWidth',5)
set(gca,'XTick',percentage_range*100)
set(gca,'XTick',percentage_range*100,'XTickLabel',{'5','','15','','25','','35','','45','','55','','65','','75','','85','','95',''})
ylim([0.9,1.1])
title('Mean square error ')
set(gca,'TickDir','out','FontSize',25,'LineWidth',2)
box off 
ylabel('mse')
xlabel('Percentage of included edges')
