%% This plots an example of the concatenated vs. averaged data 
clear all 
%% Load the data
dir='/Users/nicolehakim/Desktop/Hakim_Awh_Vogel_Rosenberg';

load([dir,'/data/compiled/oregon_raw_allsubs_averaged_vs_concat.mat']);

%% Make figure with example plot
figure; set(gcf,'Position',[57 1000 1000 1200]);

subplot(3,1,1)
plot(0:4:999,squeeze(eeg.trial_ave(1,10,:)),'LineWidth',2)
title({'Example of trial-averaged time course'})
ylabel('amplitide')
xlabel('time (ms)')
set(gca,'TickDir','out','FontSize',20,'LineWidth',1.5);
box off
ylim([-10 10])

subplot(3,1,2)
plot(0:4:999,squeeze(eeg.trial_concat(1,10,1:250)),'LineWidth',2)
title({'Example of concatenated time course','Time course from one trial'})
ylabel('amplitide')
xlabel('time (ms)')
set(gca,'TickDir','out','FontSize',20,'LineWidth',1.5);
box off
ylim([-20 20])

subplot(3,1,3)
plot(0:4:793999,squeeze(eeg.trial_concat(1,10,:)),'LineWidth',.5)
title({'Example of concatenated time course','Time course from entire experiment'})
ylabel('amplitide')
xlabel('time (ms)')
set(gca,'TickDir','out','FontSize',20,'LineWidth',1.5);
box off
ylim([-50 50])

