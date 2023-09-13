%% Isotta Rigoni
%  ~ EEG and Epilepsy Unit- Geneva HUG
% This code will calculate one connectivity matrix for each 
% one-second-epoch, for each subject


%% you will need to remove ch 31 from the rest of the data before 
%starting the preprocessing because some of Laurent's data were already average rereferenced!
clear all
close all
clc

%add fieldtrip path
addpath('C:\Users\IRIO\Desktop\Matlab toolbox\fieldtrip-master\fieldtrip-master');
addpath('func')
ft_defaults

%% variable initialisation
path='H:\Isotta\DATA\ir_mice_project\RS\';
BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\dataset';
% task='task-rest';
task='task-GS';

%define animal ID and path2data
if strcmp(task,'task-rest')
    subj=[1:33];%animal ID
    path2data=BIDSfolder;
elseif strcmp(task,'task-GS')
    subj=[12 13 15:20 22:33];%animal ID
    path2data=fullfile(BIDSfolder,'derivatives','eeg_pre_preprocessing_GS');
end

fs=4000; %sampling frequency

% prepare electrodes and neighbours for interpolation
load('H:\Isotta\DATA\ir_mice_project\RS\dataset\derivatives\elec_layout\elec_layout.mat'); %load elec layout
load('H:\Isotta\DATA\ir_mice_project\RS\dataset\derivatives\elec_layout\elec_neigh.mat'); %load elec neighbours

%% PREPROCESSING -fieldtrip-
for s=1:length(subj) %we don't have Laurent's data yet
    sub_id=subj(s);
    
    %list the sessions you have for each subj
    cnt=dir(fullfile(path2data,['sub-',sprintf('%02d',sub_id)]));
    
    %preprocess the data in each session
    for ses_idx=3:length(cnt)
        clearvars epoch_data data_tmp dataDown datFilt dataReRef dataFinal dataPreProc
        %session ID
        ses_id=cnt(ses_idx).name;
        
        raw_data_filename=fullfile(path2data,['sub-',sprintf('%02d',sub_id)],ses_id,'eeg',...
            ['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,'_eeg.mat']);
        
        if ~exist(raw_data_filename)
            continue
        else
            %load raw data
            load(raw_data_filename);
            
            %START THE PREPROCESSING------------------------
            if strcmp(task,'task-rest') %only rest data needs to be converted to FT
                %convert data to FT format
                for iep = 1:size(epoch_data,3)
                    data_tmp.trial{1,iep}=squeeze(epoch_data(:,:,iep));
                    data_tmp.time{1,iep}=0: 1/fs: size(data_tmp.trial{1,iep},2)/fs-1/fs;
                end
                data_tmp.fsample=fs;
                data_tmp.label=elec.label;
                data_tmp.dimord='chan_time';
                data_tmp.elec=elec;
                data=ft_preprocessing([], data_tmp);
            end
            if ismember(subj(s),28:33)
                %subtract channel 31 form the rest of the data -->now all data are
                %the same
                %(Guru's mice and some of Laurent's mice are already avg reref in ch 31;
                %the remaining Laurent's mice are not reref and ch 31 is set to zero)
                %it won't change Guru's data, but it will change Laurent's data
                cfg=[];
                cfg.reref='yes';
                cfg.refchannel=data.label{31, 1};
                data=ft_preprocessing(cfg,data);
            end
            % downsample to 250
            cfg.resamplefs=250;
            dataDown=ft_resampledata(cfg,data);
            
            % Butterworth filter (1-40) and notch filter
            cfg=[];
            cfg.bpfilter='yes';
            cfg.bpfreq= [1 40];
            cfg.dftfilter='yes';
            cfg.bpfiltord=6;
            dataFilt=ft_preprocessing(cfg,dataDown);
            
            % find bad channels
            bad_ch=find_bad_channels_mice(dataFilt);
            
            if ~isempty(bad_ch)
                
                % interpolate bad channels and
                cfg = [];
                cfg.badchannel=bad_ch;
                cfg.method='spline';
                cfg.neighbours=neighbours;
                data_interp = ft_channelrepair(cfg,dataFilt);
                
                % reset the original order of the channels if changed
                orig_labels=dataFilt.label;
                if ~all(strcmp(data_interp.label, orig_labels))
                    [~, idx] = ismember(orig_labels,data_interp.label);
                    data_interp.label = data_interp.label(idx);
                    data_interp.trial = cellfun(@(trl) trl(idx,:), data_interp.trial, 'uni', 0);
                end
            else
                data_interp=dataFilt;
            end
            
            %plot and check data
            %         ft_databrowser([],data_interp);title(['sub-',sprintf('%02d',sub_id),' ',ses_id])
            if isempty(bad_ch)
                msg{s,ses_idx-2}={};
            else
                msg{s,ses_idx-2}=bad_ch;
            end
            
            %average rereference
            cfg=[];
            cfg.reref='yes';
            cfg.refchannel='all';
            cfg.refmethod='avg';
            dataReRef=ft_preprocessing(cfg,data_interp);
            
            %Remove ch 31 AFTER avg reref
            cfg=[];
            cfg.channel={'all', '-31'};
            dataFinal=ft_selectdata(cfg,dataReRef);
            
            %rename final data
            dataPreProc=dataFinal;
            
            %store results in derivatives
            %define final name and folder
            final_folder=fullfile(BIDSfolder,'derivatives','eeg_preprocessing',['sub-',sprintf('%02d',sub_id)],ses_id,'eeg');
            final_filename=['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,'_eeg'];
            if ~exist(final_folder)
                mkdir(final_folder)
            end
            %save in BIDS
            save(fullfile(final_folder,final_filename),'dataPreProc');
            clear dataPreProc
        end
        %         %-------------------------- z score data and store ----------------
        %         %Convert data to matrix
        %         for ep=1:size(dataFinal.trial,2)
        %             dat(:,:,ep)=dataFinal.trial{1,ep};
        %         end
        %         %zscore the data along time
        %         dat_zscored=zscore(dat,0,2);
        %
        %         %Convert back to fieldtrip format
        %         for ep = 1:size(dataFinal.trial,2)
        %             dataFinal.trial{1,ep}=squeeze(dat_zscored(:,:,ep));
        %         end
        %
        %         %save final result
        %         dataPreProc=dataFinal;
        %
        %         %store results in derivatives
        %         %define final name and folder
        %         final_folder=fullfile(BIDSfolder,'derivatives','eeg_preprocessing_zscore',['sub-',sprintf('%02d',sub_id)],ses_id,'eeg');
        %         final_filename=['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,'_eeg'];
        %         if ~exist(final_folder)
        %             mkdir(final_folder)
        %         end
        %         %save in BIDS
        %         save(fullfile(final_folder,final_filename),'dataPreProc');
    end
end
%%
for r=1:size(msg,1)
    for c=1:size(msg,2)
        if ~isempty(msg{r,c})
            disp(['in M',num2str(subj(r)),' ch ', char(msg{r,c}),' was removed'])
            
        end
    end
end

