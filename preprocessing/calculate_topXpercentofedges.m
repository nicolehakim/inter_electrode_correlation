%% This plots the top ten percent of edges 
clear all; clc
%% Load the data :)
load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/BigStudy/Data/fcmat_allsubs_RI_and_baseline_filtered.mat');%Oregon
chicago=load('/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/fc_mat_uchicagosubs_n165.mat');%Chicago

%% Settings 
percentage_range=[.02,.06,.1,.14,.15,.18,.2,.22,.26,.3,.34,.38,.42,.46,.5];%the percentage of electrodes that we will include in this analysis

nElec=length(bs.fc_mat.overlappingElecs.elecs);

nEdges=((nElec*nElec)-nElec)/2;

for pp=1:length(percentage_range)
    percentage=percentage_range(pp);
    % Get the to X%
    nElec2include=ceil(nEdges*(percentage));
    if mod(nElec2include,2)~=0%if there's not an even number of edges, then make it even
        nElec2include=nElec2include+1;
    end
    
    rhoEdges.color=nan(nElec,nElec);
    rhoEdges.shape=nan(nElec,nElec);
    for ed1=1:nElec
        for ed2=1:nElec
            %train color
            [rhoEdges.color(ed1,ed2),p]=corr(squeeze(bs.fc_mat.overlappingElecs.RI.color(:,ed1,ed2)),bs.fc_mat.k.color);
            
            %train shape
            [rhoEdges.shape(ed1,ed2),p]=corr(squeeze(bs.fc_mat.overlappingElecs.RI.shape(:,ed1,ed2)),bs.fc_mat.k.shape);
        end
    end
    
    sigEdges.color=zeros(nElec*nElec,1);
    sigEdges.shape=zeros(nElec*nElec,1);
    
    %color
    tmp=sort(rhoEdges.color(~isnan(rhoEdges.color)),1,'ascend');
    sigEdges.color(ismember(rhoEdges.color,tmp(1:nElec2include)))=-1;%negative edges
    sigEdges.color(ismember(rhoEdges.color,tmp(end-nElec2include+1:end)))=1;%positive edges
    sigEdges.color=reshape(sigEdges.color,[nElec,nElec]);
    
    %shape
    tmp=sort(rhoEdges.shape(~isnan(rhoEdges.shape)),1,'ascend');
    sigEdges.shape(ismember(rhoEdges.shape,tmp(1:nElec2include)))=-1;%negative edges
    sigEdges.shape(ismember(rhoEdges.shape,tmp(end-nElec2include+1:end)))=1;%positive edges
    sigEdges.shape=reshape(sigEdges.shape,[nElec,nElec]);
    
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1<=ed2
                sigEdges.shape(ed1,ed2)=0;
                sigEdges.color(ed1,ed2)=0;
            end
        end
    end
    
    %% Plot Oregon: color & shape separately
    figure; set(gcf,'Position',[57 1000 800 1400]);
    
    %color
    subplot(2,1,1)
    imagesc(sigEdges.color)
    colormap([76 160 174; 255 255 255; 240 157 56]/255)
    title({['Top ',num2str(percentage*100),'% of edges'],'Oregon: color'})
    set(gca,'YTick',[1:nElec],'YTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    set(gca,'XTick',[1:nElec],'XTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    xtickangle(45)
    
    %shape
    subplot(2,1,2)
    imagesc(sigEdges.shape)
    colormap([76 160 174; 255 255 255; 240 157 56]/255)
    title('Oregon: shape')
    set(gca,'YTick',[1:nElec],'YTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    set(gca,'XTick',[1:nElec],'XTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    xtickangle(45)
    
    %% Get Chicago top x% of edges
    rhoEdges_chicago=nan(nElec,nElec);
    for ed1=1:nElec
        for ed2=1:nElec
            %train color
            [rhoEdges_chicago(ed1,ed2),p]=corr(squeeze(chicago.fc_mat.eeg(:,ed1,ed2)),chicago.fc_mat.beh);
        end
    end
    
    sigEdges.chicago=zeros(nElec*nElec,1);
    
    tmp=sort(rhoEdges_chicago(~isnan(rhoEdges_chicago)),1,'ascend');
    sigEdges.chicago(ismember(rhoEdges_chicago,tmp(1:nElec2include)))=-1;%negative edges
    sigEdges.chicago(ismember(rhoEdges_chicago,tmp(end-nElec2include+1:end)))=1;%positive edges
    sigEdges.chicago=reshape(sigEdges.chicago,[nElec,nElec]);
    
    for ed1=1:nElec
        for ed2=1:nElec
            if ed1<=ed2
                sigEdges.chicago(ed1,ed2)=0;
            end
        end
    end
    
    %% Plot the significant edges for Oregon & Chicago & their overlap
    figure; set(gcf,'Position',[57 1000 2400 600]);
    %overlap of Oregon color & shapep
    subplot(1,3,1)
    sigEdges_overlap.color_shape=zeros(nElec,nElec);
    sigEdges_overlap.color_shape(sigEdges.color==1 & sigEdges.shape==1)=1;
    sigEdges_overlap.color_shape(sigEdges.color==-1 & sigEdges.shape==-1)=-1;
    imagesc(sigEdges_overlap.color_shape)
    colormap([76 160 174; 255 255 255; 240 157 56]/255)
    set(gca,'YTick',[1:nElec],'YTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    set(gca,'XTick',[1:nElec],'XTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    xtickangle(45)
    title({'Oregon'})
    
    %chicago
    subplot(1,3,2)
    imagesc(sigEdges.chicago)
    colormap([76 160 174; 255 255 255; 240 157 56]/255)
    set(gca,'YTick',[1:nElec],'YTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    set(gca,'XTick',[1:nElec],'XTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    xtickangle(45)
    title({['Top ',num2str(percentage*100),'% of edges'],'Chicago'})
    
    %chicago oregon overlap
    subplot(1,3,3)
    sigEdges_overlap.oregon_chicago=zeros(nElec,nElec);
    sigEdges_overlap.oregon_chicago(sigEdges_overlap.color_shape==1 & sigEdges.chicago==1)=1;
    sigEdges_overlap.oregon_chicago(sigEdges_overlap.color_shape==-1 & sigEdges.chicago==-1)=-1;
    imagesc(sigEdges_overlap.oregon_chicago)
    colormap([76 160 174; 255 255 255; 240 157 56]/255)
    set(gca,'YTick',[1:nElec],'YTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    set(gca,'XTick',[1:nElec],'XTickLabel',bs.fc_mat.overlappingElecs.elecs,'fontsize', 25)
    xtickangle(45)
    title({'Oregon/Chicago overlap'})
    
    %% Save edges
    save(['/Users/nicolehakim/Desktop/Maude/Experiments/Experiments_MRI/UChicago_studies_compiled/data/oregon_chicago_top',num2str(percentage*100),'percedges.mat'],'sigEdges_overlap','sigEdges')
    
end

