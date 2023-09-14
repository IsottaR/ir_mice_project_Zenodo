%% Isotta Rigoni
%  ~ EEG and Epilepsy Unit- Geneva HUG
%AVERAGE CONNECTOMES BETWEEN DAY 28 AND 29

clear all
close all
clc

%add path
addpath('C:\Users\IRIO\Desktop\PROJECTS\NeDepi\ieds\graph_signal_processing_functions')
% ft_defaults

%% variable initialisation
%-----------------------------------------
day={'d28' 'd29'}; %time points of observation
subj=[12 13 15:20 22:33]; %animal ID
% %-----------------------------------------

BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\data2publish\derivatives';
task='task-rest';
% n_elec=30;
derivative_folder='wpli';

%% calculate

for s=1:length(subj)
    %subject id
        sub_id=subj(s);
    for d=1:length(day)
        %get network metric for each session: d0, d28 or d29
        
        %session ID
        ses_id=['ses-',char(day(d))];
        
        % ------------ Load network data ----------------------------------

        %load connectivity matrix
        load(fullfile(BIDSfolder,derivative_folder,['sub-',sprintf('%02d',sub_id)],ses_id,'eeg',...
            ['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,'.mat']));
        
        %save metrics
        wpli_avg_delta_2d(:,:,d)=wpli_avg_delta;
        wpli_avg_lowTheta_2d(:,:,d)=wpli_avg_lowTheta;
        wpli_avg_highTheta_2d(:,:,d)=wpli_avg_highTheta;
        wpli_avg_beta_2d(:,:,d)=wpli_avg_beta;
        wpli_avg_gamma_2d(:,:,d)=wpli_avg_gamma;
        wpli_avg_2d(:,:,d)=wpli_avg;
        
        clear wpli_avg_delta wpli_avg_lowTheta wpli_avg_highTheta wpli_avg_beta wpli_avg_gamma
    end
    
    wpli_avg_delta=squeeze(mean(wpli_avg_delta_2d,3));
    wpli_avg_lowTheta=squeeze(mean(wpli_avg_lowTheta_2d,3));
    wpli_avg_highTheta=squeeze(mean(wpli_avg_highTheta_2d,3));
    wpli_avg_beta=squeeze(mean(wpli_avg_beta_2d,3));
    wpli_avg_gamma=squeeze(mean(wpli_avg_gamma_2d,3));
    wpli_avg=squeeze(mean(wpli_avg_2d,3));
    
    %store results in derivatives
    %define final name and folder
    final_folder=fullfile(BIDSfolder,'wpli',['sub-',sprintf('%02d',sub_id)],'ses-d28-d29','eeg');
    final_filename=['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,''];
    if ~exist(final_folder)
        mkdir(final_folder)
    end
    %save in BIDS
    save(fullfile(final_folder,final_filename),...
        'wpli_avg','wpli_avg_delta','wpli_avg_lowTheta','wpli_avg_highTheta','wpli_avg_beta','wpli_avg_gamma','-v7.3')
    
    clear wpli_avg_delta wpli_avg_lowTheta wpli_avg_highTheta wpli_avg_beta wpli_avg_gamma
    
end