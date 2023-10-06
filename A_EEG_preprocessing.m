%% Isotta Rigoni
%  ~ EEG and Epilepsy Unit- Geneva HUG
% This code will import epochs of raw Eepicranial EEG and preprocess them
%The preprocessing is performed Fieldtrip and is described in the paragraph
%"EEG recordings and analyses"

% -------> change path at line 14 and 19

clear all
close all
clc

%add fieldtrip path
addpath('C:\Users\IRIO\Desktop\Matlab toolbox\fieldtrip-master\fieldtrip-master');
addpath('func')
ft_defaults

%% variable initialisation
BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\data2publish'; %insert your path to the data
task='task-rest';

%list of animals 2 analyse
cnt=dir(BIDSfolder);

fs=4000; %sampling frequency

% load electrodes and neighbours for interpolation
load(fullfile(BIDSfolder,'derivatives\elec_layout\elec_layout31.mat')); %load elec layout
load(fullfile(BIDSfolder,'derivatives\elec_layout\elec_neigh.mat')); %load elec neighbours

%% PREPROCESSING -fieldtrip-
for s=4:size(cnt,1) 
    sub_id=cnt(s).name;
    
    %list the sessions you have for each subj
    cnt_ses=dir(fullfile(BIDSfolder,sub_id));
    
    %preprocess the data in each session
    for ses_idx=3:length(cnt_ses)
        clearvars epoch_data data_tmp dataDown datFilt dataReRef dataFinal dataPreProc
        %session ID
        ses_id=cnt_ses(ses_idx).name;
        
        raw_data_filename=fullfile(BIDSfolder,sub_id,ses_id,'eeg',...
            [sub_id,'_',ses_id,'_',task,'_eeg.mat']);
        
        if ~exist(raw_data_filename)
            continue
        else
            %load raw data
            load(raw_data_filename);
            
            %START THE PREPROCESSING------------------------
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
            
            if ismember(str2num(sub_id(5:end)),28:33)
                %the raw data of these animals were avg re-referenced to 
                %channel 31, so subtract channel 31-->now all data are the same
                
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
            final_folder=fullfile(BIDSfolder,'derivatives','eeg_preprocessing',sub_id,ses_id,'eeg');
            final_filename=[sub_id,'_',ses_id,'_',task,'_eeg'];
            if ~exist(final_folder)
                mkdir(final_folder)
            end
            %save in BIDS
            save(fullfile(final_folder,final_filename),'dataPreProc');
            clear dataPreProc
        end
    end
end

%% check if any channel was removed from any animal
for r=1:size(msg,1)
    for c=1:size(msg,2)
        if ~isempty(msg{r,c})
            disp(['in M',num2str(subj(r)),' ch ', char(msg{r,c}),' was removed'])
            
        end
    end
end

