%% Isotta Rigoni
%  ~ EEG and Epilepsy Unit- Geneva HUG

%With this script, you can run the statistical analyses between:
% - graph measures at day 0 vs day 28 (group A)
% - graph measures at day 28 vs day 29(group B)

%The analyses are run for GLOBAL (GE, GCC and LI), 
%HEMISPHERIC (HCC and HE, ipsilateral and contralateral) 
% and NODAL measures (NE and CC). The LI is also compared 
%with a zero-mean distribution (line 203)

% -------> change path at line 21 and 33

clear all
close all
clc

%add path
addpath('func')
addpath('C:\Users\IRIO\Desktop\Matlab toolbox\fieldtrip-master\fieldtrip-master')
ft_defaults

%% variable initialisation
%decide which group to compare
day={'d0' 'd28'}; %time points of observation
subj=[1:7 12 14:17 21:23 ]; %animal ID
%-----------------------------------------
% day={'d28' 'd29'}; %time points of observation
% subj=[12 13 15:20 22:33]; %animal ID
%-----------------------------------------

BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\data2publish\derivatives';
task='task-rest';
n_elec=30;
var_labels_glob_hemisp={'GE','LI','avgCC','GE_L','GE_R','CC_R','CC_L'};%global and hemispheric measures
var_labels_nodal={'CC','NE'};%nodal measures

derivative_folder='wpli';

%---------load electrode layout for label 
load(fullfile(BIDSfolder, 'elec_layout', 'elec_layout30.mat'))
%% LOAD and ORGANISE the data
for d=1:length(day)
    %get network metric for each session: d0, d28 or d29
    %session ID
    ses_id=['ses-',char(day(d))];
    for s=1:length(subj)
        clearvars GE STR LI CC
        
        %subject id
        sub_id=subj(s);
        
        % ------------ Load network data ----------------------------------
        load(fullfile(BIDSfolder,['network_metrics_',derivative_folder],['sub-',sprintf('%02d',sub_id)],ses_id,'eeg',...
            ['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,'_network_metrics.mat']));
        
        %-----------global measures
        GE_2_stat(d,s,:)=squeeze(GE);
        LI_2_stat(d,s,:)=squeeze(LI);
        avgCC_2_stat(d,s,:)=avgCC;
        
        %-----------hemispheric measures
        
        GE_R_2_stat(d,s,:)=squeeze(GE_hemi(1,:));
        GE_L_2_stat(d,s,:)=squeeze(GE_hemi(2,:));
        
        CC_R_2_stat(d,s,:)=squeeze(CC_hemi(1,:));
        CC_L_2_stat(d,s,:)=squeeze(CC_hemi(2,:));
        
        %-----------single node measures
        %(store them in Fieldtrip-like PSD structure to run cluster-based permutation test)
        eval(['CC_2_stat',char(day(d)),'.powspctrm(s,:,:)=squeeze(CC);']);
        eval(['NE_2_stat',char(day(d)),'.powspctrm(s,:,:)=squeeze(NE);']);
        
    end
    %complete CC and NE structures as if they were PSD 
    str='subj_chan_freq';
    eval(['CC_2_stat',char(day(d)),'.label=arrayfun(@num2str, [1:30]'', ''UniformOutput'', 0);']);
    eval(['CC_2_stat',char(day(d)),'.freq=[1 2 3 4 5 6];']);
    eval(['CC_2_stat',char(day(d)),'.dimord=str;']);
    
    eval(['NE_2_stat',char(day(d)),'.label=arrayfun(@num2str, [1:30]'', ''UniformOutput'', 0);']);
    eval(['NE_2_stat',char(day(d)),'.freq=[1 2 3 4 5 6];']);
    eval(['NE_2_stat',char(day(d)),'.dimord=str;']);
end

 band={'delta','lowTheta','highTheta','beta','gamma',''}

%% COMPARISONS
%% -----------------compare GLOBAL and HEMISPHERIC measures-----------------------
%GE, LI, GE_L, GE_R avgCC CC_L CC_R
%--> make Fig 3
for v=1:length(var_labels_glob_hemisp)
    eval(['var=',char(var_labels_glob_hemisp(v)),'_2_stat;']);
    for b=1:5%size(band,2)
        p(b)=signrank(squeeze(var(2,:,b)), squeeze(var(1,:,b)));
    end
    
    %Bonferroni correction over 5 freq bands
    p=p*5;
    for b=1:5%size(band,2)
        if p(b)<.05
            diff=median(var(2,:,b))-median(var(1,:,b));
            disp(['freq: ',char(band(b)),'; ',char(var_labels_glob_hemisp(v)),'(',char(day(2)),') - ',char(var_labels_glob_hemisp(v)),'(',char(day(1)),') = ',num2str(diff),'; p = ',num2str(p(b))])
            %plot
            fig=scatter_box_plot(squeeze(var(1,:,b)),squeeze(var(2,:,b)),day);
            title([char(var_labels_glob_hemisp(v)),'(',char(day(2)),') - ',char(var_labels_glob_hemisp(v)),'(',char(day(1)),') = ',num2str(diff),'; p = ',num2str(p(b))])
            xlabel('time point')
            ylabel([char(var_labels_glob_hemisp(v)), ' in ',char(band(b))])
        end
    end
end

%% -----------------compare NODAL measures with cluster-based permutest-----------------------------------
%CC and NE --> only look in lowTheta and highTheta because there are no 
%significant results in the other bands
%--> make Fig 4
band2analyze=[2 3];

%---------prepare neighbours
cfg_neigh = [];
cfg_neigh.layout = layout;
cfg_neigh.method = 'triangulation';
cfg_neigh.compress = 'no';
cfg_neigh.feedback = 'yes';
neighbours = ft_prepare_neighbours(cfg_neigh);
clear layout

%load layout to PLOT (with mice head)
load(fullfile(BIDSfolder, 'elec_layout','mouse_layout_modif.mat'))

for v=1:length(var_labels_nodal)
    clear p
    
    eval(['var1=',char(var_labels_nodal(v)),'_2_stat',char(day(1))]);
    eval(['var2=',char(var_labels_nodal(v)),'_2_stat',char(day(2))]);
    n_sub=size(var1.powspctrm,1);%nsub is the same ni both variable, it's a WHITHIN-SUB
    
    %calculate average difference between the two powspctrm
    cfg=[];cfg.operation='subtract';cfg.parameter='powspctrm';
    var=ft_math(cfg,var2,var1);
    
    
    for b=1:size(band2analyze,2)
        
        clear stat2
        %prepare cfg for cluster based permutation test
        cfg                 = [];
        cfg.channel          = 'all'; cfg.avgovergchan = 'no';
        cfg.frequency        = [band2analyze(b) band2analyze(b)]; cfg.avgoverfreq = 'yes';
        cfg.parameter        = 'powspctrm';
        cfg.method           = 'ft_statistics_montecarlo'; cfg.statistic='ft_statfun_depsamplesT';
        cfg.correctm         = 'cluster'; cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum';
        cfg.clusterthreshold = 'nonparametric_common';
        cfg.minnbchan        = 2;
        cfg.tail             = 0; cfg.clustertail= cfg.tail; %two sides test
        cfg.alpha            = 0.05; cfg.correcttail = 'alpha';cfg.computeprob = 'yes'; cfg.numrandomization = 5000;
        cfg.neighbours       = neighbours;
        
        %design matrix
        design = zeros(2, n_sub*2);
        design(1,:) = [1:n_sub 1:n_sub];
        design(2,:) = [ones(1,n_sub) ones(1,n_sub)*2];
        cfg.design = design;
        cfg.uvar   = 1;
        cfg.ivar   = 2;
        
        %do actual test
        stat2 = ft_freqstatistics(cfg, var1, var2);
        
        %find the pvalues in stat2.posclusters or stat2.negclusters, and
        %correct it by 4
        
        %plot if any channel is significant
%         if any(stat2.mask)
            
            %------------ PLOT TOPOPLOT -----------------------
            cfg = [];
            cfg.parameter= 'powspctrm';
            cfg.xlim=[band2analyze(b) band2analyze(b)];
            cfg.channel=var.label;
            cfg.marker='on';
            cfg.markersize=3;
            cfg.layout= lay_modif;
            cfg.highlight='on';
            cfg.highlightsymbol='o';
            cfg.highlightcolor=[1 1 0];
            cfg.highlightchannel=stat2.label(stat2.mask(:,1));%plot the ch that are signif
            cfg.zlim='maxabs';
            
            fig=figure;
            ft_topoplotER(cfg, var);
            
            title([char(var_labels_nodal(v)),' in ',char(band(band2analyze(b)))]);
            ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
            colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
            colorbar
            clear fig

%         end
    end
end

%% compare LI with zero at d0 and d28 in theta and alpha
var_labels3={'LI'}
for b=2:3
    for v=1:length(var_labels3)
        eval(['var=',char(var_labels3(v)),'_2_stat;']);
        for d=1:length(day)
            p=signrank(squeeze(var(d,:,b)));
            
            p=p*4;%Bonferroni correction for (2 days)*(2 freq bands)
            if p<.05
                disp(['freq: ',char(band(b)),'; ',char(var_labels3(v)),'; p = ',num2str(p)])
                %plot
                fig=figure
                boxplot(squeeze(var(d,:,b)));
                title([char(var_labels3(v)),'in ',char(band(b)),'; p = ',num2str(p)])
                xlabel(char(day(d)))
                ylabel([char(var_labels3(v)), ' in ',char(band(b))])
                clear fig
            end
            clear p
        end
    end
end

