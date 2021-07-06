%% Check that the data in the two experiments are correctly aligned
clear all 
%% Load the big study data 
dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg/';
bs_filename=[dir,'/data/compiled/','oregon_site_fc_mat.mat'];
bs=load(bs_filename);

%% Load psych science data 
ps=load([dir,'data/compiled/chicago_site_fc_mat.mat']);

%% Find where they overlap
bs_e_idx=nan(1,length(ps.fc_mat.elecOrder));
for elec=1:length(ps.fc_mat.elecOrder)
    bs_e_idx(elec)=find(contains(bs.fc_mat.elecs_new,ps.fc_mat.elecOrder(elec)));
end

save([dir,'/data/preprocessed/oregon/Oregon_electrode_order.mat'],'bs_e_idx')
%% Order the big study data, so that it matches the psych science data order 
bs.fc_mat.overlappingElecs.baseline.color=bs.fc_mat.baseline.color(:,bs_e_idx,bs_e_idx);
bs.fc_mat.overlappingElecs.baseline.shape=bs.fc_mat.baseline.shape(:,bs_e_idx,bs_e_idx);
bs.fc_mat.overlappingElecs.baseline.bothTasks=bs.fc_mat.baseline.bothTasks(:,bs_e_idx,bs_e_idx);

bs.fc_mat.overlappingElecs.RI.color=bs.fc_mat.RI.color(:,bs_e_idx,bs_e_idx);
bs.fc_mat.overlappingElecs.RI.shape=bs.fc_mat.RI.shape(:,bs_e_idx,bs_e_idx);
bs.fc_mat.overlappingElecs.RI.bothTasks=bs.fc_mat.RI.bothTasks(:,bs_e_idx,bs_e_idx);
%% save the names of the overlapping electrodes to the big study file also 
bs.fc_mat.overlappingElecs.elecs=ps.fc_mat.elecOrder;
%% Save the new reconfigured data!
save(bs_filename,'bs','-v7.3');

