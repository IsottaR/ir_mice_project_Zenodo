%% Isotta Rigoni
%  ~ EEG and Epilepsy Unit- Geneva HUG

%This code calculates WPLI connectivity on epicranial eeg data
%For each animal and each session (d0, d28 and d29), we calculate one WPLI
%matrix and then split it into frequency bands, as described in the
%paragraph "Connectivity analyses"

% -------> change path at line 16 and 21

clear all
close all
clc

%add fieldtrip path
addpath('C:\Users\IRIO\Desktop\Matlab toolbox\fieldtrip-master\fieldtrip-master');

ft_defaults

%% variable initialisation
BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\data2publish';
task='task-rest';

%list of animals 2 analyse
cnt=dir(fullfile(BIDSfolder,'derivatives','eeg_preprocessing'));

freqs=linspace(1,40,100);
freq_band_labels={'delta','theta','alpha','beta','gamma','broadband'};
band_freq=[1 4; 4 8; 8 12; 12 30; 30 40; 1 40];
%% PREPROCESSING -fieldtrip-
for s=3:size(cnt,1)
    sub_id=cnt(s).name;
    
    %list the sessions you have for each subj
    cnt_ses=dir(fullfile(BIDSfolder,'derivatives','eeg_preprocessing',sub_id));
    
    %preprocess the data in each session
    for ses_idx=3:length(cnt_ses)
        clearvars freq  trials results wpli wpli_avg_delta wpli_avg_theta wpli_avg_alpha wpli_avg_beta wpli_avg_gamma
        
        %session ID
        ses_id=cnt_ses(ses_idx).name;
        
        %name of the preprocessed data file
        raw_data_filename=fullfile(BIDSfolder,'derivatives','eeg_preprocessing',sub_id,ses_id,'eeg',...
            [sub_id,'_',ses_id,'_',task,'_eeg.mat']);
        
        if ~exist(raw_data_filename)
            continue
        else
            
            
            %load data, if it exists
            load(raw_data_filename)
            
            %get fourier spectrum
            cfg             = [];
            cfg.output      = 'fourier';
            cfg.taper       = 'dpss'; % or dpss
            cfg.method      = 'mtmfft';
            cfg.foi         = freqs;
            cfg.pad         = 'nextpow2';
            cfg.keeptrials  = 'yes';
            if strcmp(cfg.taper, 'dpss')
                cfg.tapsmofrq   = 1; % number of tapers
            end
            freq     = ft_freqanalysis(cfg, dataPreProc);
            
            % debiased weighted phase lag index
            cfg         	= [];
            cfg.method      = 'wpli_debiased';
            cfg.feedback = 'none';
            results = ft_connectivityanalysis(cfg, freq);
            wpli = results.([cfg.method 'spctrm']).*~eye(size(dataPreProc.trial{1, 1},1));
            %put negative values to zero
            wpli = wpli.*(wpli > 0);
            %put NaN to zero
            wpli(isnan(wpli))=0;
            
            
            %get PDC matrixes in different frequency bands
            wpli_avg=squeeze(mean(wpli,3)); % averaging in frequency [d x d]
            
            f1_delta= find(results.freq  >=1,1,'first');
            f2_delta= find(results.freq  >=4,1,'first');
            wpli_avg_delta=squeeze(mean(wpli(:,:,f1_delta:f2_delta),3)); % averaging in frequency [d x d]
            
            f1_theta= find(results.freq  >=4,1,'first');
            f2_theta= find(results.freq  >=8,1,'first');
            wpli_avg_lowTheta=squeeze(mean(wpli(:,:,f1_theta:f2_theta),3)); % averaging in frequency [d x d]
            
            f1_alpha= find(results.freq  >8,1,'first');
            f2_alpha= find(results.freq  >=12,1,'first');
            wpli_avg_highTheta=squeeze(mean(wpli(:,:,f1_alpha:f2_alpha),3)); % averaging in frequency [d x d ]
            
            f1_beta= find(results.freq  >=12,1,'first');
            f2_beta= find(results.freq  >=30,1,'first');
            wpli_avg_beta=squeeze(mean(wpli(:,:,f1_beta:f2_beta),3)); % averaging in frequency [d x d ]
            
            f1_gamma= find(results.freq  >=30,1,'first');
            f2_gamma= find(results.freq  >=40,1,'first');
            wpli_avg_gamma=squeeze(mean(wpli(:,:,f1_gamma:f2_gamma),3)); % averaging in frequency [d x d ]
            
            %store results in derivatives
            %define final name and folder
            final_folder=fullfile(BIDSfolder,'derivatives','wpli',sub_id,ses_id,'eeg');
            final_filename=[sub_id,'_',ses_id,'_',task,''];
            if ~exist(final_folder)
                mkdir(final_folder)
            end
            %save in BIDS
            save(fullfile(final_folder,final_filename),...
                'wpli_avg','wpli_avg_delta','wpli_avg_lowTheta','wpli_avg_highTheta','wpli_avg_beta','wpli_avg_gamma','-v7.3')
        end
        
    end
end
