%% Isotta Rigoni
%  ~ EEG and Epilepsy Unit- Geneva HUG

clear all
close all
clc

%add BCT path
addpath('C:\Users\IRIO\Desktop\Matlab toolbox\BCT\2019_03_03_BCT')

%% variable initialisation
BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\data2publish';
task='task-rest';
derivative_folder='wpli';

%list of animals 2 analyse
cnt=dir(fullfile(BIDSfolder,'derivatives','eeg_preprocessing'));

% electrodes
elec_2remove=[30 28 29 26 27 25 24 23 21 19 17 20 18 16 14 7; 15 12 13 10 11 22 8 9 3 5 6 1 2 4 14 7];
elec_2keep=[1 2 3 4 5 6 8 9 10 11 12 13 15 22; 16 17 18 19 20 21 23 24 25 26 27 28 29 30]; %remove elec 7 and 14 form both hemi

elec_2keep_bool=zeros(2,30);
elec_2keep_bool(1,elec_2keep(1,:))=1;
elec_2keep_bool(2,elec_2keep(2,:))=1;
elec_2keep_bool=logical(elec_2keep_bool)
hemisp_label={'R','L'};

% Selection of frequency band
band= [{'_delta'}, {'_lowTheta'},{'_highTheta'} ,{'_beta'} ,{'_gamma'}, {''}];

%% NETWORK ANALYSES
for s=3:size(cnt,1)
    sub_id=cnt(s).name;
    
    %list the sessions you have for each subj
    cnt_ses=dir(fullfile(BIDSfolder,'derivatives','eeg_preprocessing',sub_id));
    
    %get network metric for each session
    for ses_idx=3:length(cnt_ses)
        clearvars GE SO SO_h LI PDC BC
        %session ID
        ses_id=cnt_ses(ses_idx).name;
        
        % Load WPLI data
        load(fullfile(BIDSfolder,'derivatives',derivative_folder,sub_id,ses_id,'eeg',...
            [sub_id,'_',ses_id,'_',task,'.mat']));
        
        for iband = 1:size(band,2)
            eval(['Data_matrix = wpli_avg' band{iband} ';'])
            
            %% ------ NODAL MEASURES
            
            %nodal efficiency (inverse of shortest path of each node)
            A=Data_matrix;
            A_inv=1./(A);%invert the matrix --> get length
            A_inv(isinf(A_inv))=0;%put inf elem to zero
            NE(:,iband)=sum(distance_inv_wei(A_inv),2)./(size(A_inv,1) - 1);
            
            %clustering coefficient of each node
            CC(:,iband) = clustering_coef_wu(Data_matrix);
            
            %STRENGTH
            STR(:,iband)  = strengths_und(Data_matrix);
            
            %% ------- GLOBAL MEASURES
            %global efficiency of the network
            GE(iband) = efficiency_wei(Data_matrix,0);
            
            avgCC(iband)=mean(CC(:,iband));
            
            %% ------- HEMISPHERIC MEASURES
            %put sensor 7 and 14 to zero
            Data_matrix=Data_matrix*diag([1 1 1 1 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
            Data_matrix=(Data_matrix'*diag([1 1 1 1 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]));
            
            %get total strength of each node (not counting 7 and 14)
            total_STR=strengths_und(Data_matrix)';
            
            %get GE for each HEMISPHERE
            for h=1:length(hemisp_label)
                clearvars Data_hemi
                
                %keep electrodes of current hemisphere
                Data_hemi=Data_matrix(elec_2keep_bool(h,:),elec_2keep_bool(h,:));
                
                %global efficiency of the network
                GE_hemi(h,iband) = efficiency_wei(Data_hemi,0); %Efficiency only considering intrahemisph connections
                
                %clutering coefficient of the network
                CC_hemi(h,iband) = mean(clustering_coef_wu(Data_hemi));%Clustering only considering intrahemisph connections
                
                %define electrodes 2 keep
                electrode_labels{1,h}=find(elec_2keep_bool(h,:));
                
                %strength of the nodes in the current hemisph
                STR_h(h,iband)=sum(STR(elec_2keep_bool(h,:),iband));
                
                clear Data_hemi
            end
            
            %LI (laterality index)--> counting ALL connections (even
            %those with elec 7 and 14), but elec 7 and 14 are removed
            %from the hemispheres
            sumLH=STR_h(2,iband);
            sumRH=STR_h(1,iband);
            sumALL=sum(STR(:,iband),1);
            LI(1,iband)=(sumLH-sumRH)./sumALL;
            
            clear Data_matrix
        end
        
        %store results in derivatives
        %define final name and folder
        final_folder=fullfile(BIDSfolder,'derivatives',['network_metrics_',derivative_folder],sub_id,ses_id,'eeg');
        final_filename=[sub_id,'_',ses_id,'_',task,'_network_metrics'];
        if ~exist(final_folder)
            mkdir(final_folder)
        end
        %save in BIDS
        save(fullfile(final_folder,final_filename),'GE','STR','avgCC','LI','CC_hemi','GE_hemi','NE','CC','band','hemisp_label','electrode_labels','-v7.3');
        
    end
end


