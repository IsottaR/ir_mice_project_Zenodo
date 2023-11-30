%% Isotta Rigoni
%  ~ EEG and Epilepsy Unit- Geneva HUG

%This script investigate the level of variability in adjecency matrices (all frequency
%bands) between day 28 and day 29 and in epileptic activities (EA). It runs
%the analyses described in the second paragraph on the section "Statistical analyses"

% -------> change path at line 21

clear all
close all
clc

addpath('func')
%% subject initialisation
day={'d28' 'd29'}; %time points of observation
subj=[12 13 15:20 22:33]; %animal ID
%--------------------------------------

%% variable initialisation
BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\data2publish';
task='task-rest';
derivative_folder='wpli';
f_band_name=[ {'delta'}, {'lowTheta'},{'highTheta'} ,{'beta'} ,{'gamma'}, {'broadband'}];

%% load EA info
EAinfo28=readtable(fullfile(BIDSfolder,'EA_info.xlsx'),'Sheet','d28'); %change path
EAinfo29=readtable(fullfile(BIDSfolder,'EA_info.xlsx'),'Sheet','d29'); %change path
var_labels3=EAinfo28.Properties.VariableNames(2:8);

%% load and organize epileptic activities (EA) markers
for s=1:length(subj)

%       subject id
        sub_id=subj(s);

        %------------ get epilepsy measures of each animal --------------
        sub_idx28=find(table2array(EAinfo28(:,1))==sub_id);
        measures2plot(1,s,:)=table2array(EAinfo28(sub_idx28,2:8));
        
        sub_idx29=find(table2array(EAinfo29(:,1))==sub_id);
        measures2plot(2,s,:)=table2array(EAinfo29(sub_idx29,2:8));
        
        subjects28(s,1)=table2array(EAinfo28(sub_idx28,1));
        subjects29(s,1)=table2array(EAinfo29(sub_idx29,1));

        if (subjects28(s,1)~=subj(1,s)) || (subjects29(s,1)~=subj(1,s))
            error('the subjects are not consistent')
        end
end

%% load and organise connectomes 
for d=1:length(day)
    %get network metric for each session: d0, d28 or d29
    %session ID
    ses_id=['ses-',char(day(d))];
    for s=1:length(subj)
        clearvars GE STR LI CC
        
        %subject id
        sub_id=subj(s);
        
        %load connectivity matrix
        load(fullfile(BIDSfolder,'derivatives',derivative_folder,['sub-',sprintf('%02d',sub_id)],ses_id,'eeg',...
            ['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,'.mat']));
        
        %save metrics
        wpli_all(:,:,:,d,s)=cat(3,wpli_avg_delta, wpli_avg_lowTheta,...
            wpli_avg_highTheta,wpli_avg_beta,wpli_avg_gamma,wpli_avg);
    end
end


%% do ICC analyses on CONNECTION VALUES (functional connectomes) -> make Fig 5.a
for b=1:size(wpli_all,3)-1
    %average across animals
    clear curr_mat_x curr_mat_y var 1 var2

    curr_mat_x=mean(squeeze(wpli_all(:,:,b,1,:)),3);%day28
    curr_mat_y=mean(squeeze(wpli_all(:,:,b,2,:)),3);%day29
    
    %vectorize upper matrix
    var1=curr_mat_x(triu(true(size(curr_mat_x)),1));
    var2=curr_mat_y(triu(true(size(curr_mat_y)),1));
        
    %calculate the intra-class correlation coefficient
    %--> ICC for AVERAGE measures
    [r(b), LB(b), UB(b), F(b), df1(b), df2(b), p(b)] = ICC([var1 var2], 'A-k');
    
    report_string=['ICC=',num2str(r(b)),'[',num2str(LB(b)),'-',num2str(UB(b)),'],F=',num2str(F(b)),'(',num2str(df1(b)),...
        ',',num2str(df2(b)),'), p=',num2str(p(b)*(size(wpli_all,3)-1)) ];
    
    %display results
    fig=figure;
    scatter(var1,var2)
    xlabel(['avg ',char(f_band_name(b)),' FC, ',char(day(1))])
    ylabel(['avg ',char(f_band_name(b)),' FC, ',char(day(2))])
    title(report_string,'Interpreter','none') %Bonferroni-correct
    disp(['avg ',char(f_band_name(b)),' FC, ', report_string]);
    
end

clearvars r LB UB F df1 df2 p
%% do correlation analyses on epileptic activities EA (saeizures etc.) -> make rest of Fig 5

for v=1:length(var_labels3)
    %---------------------- (metric28+metric29)/2 ----------------------
    eval(['var1=squeeze(measures2plot(1,:,v));']);%day 28
    eval(['var2=squeeze(measures2plot(2,:,v));']);%day 29
        
     %calculate the intra-class correlation coefficient
    %--> ICC for SINGLE measures
    [r, LB, UB, F, df1, df2, p] = ICC([var1; var2]', 'A-1');
    report_string=['ICC=',num2str(r),'[',num2str(LB),'-',num2str(UB),'],F=',num2str(F),'(',num2str(df1),...
        ',',num2str(df2),'), p=',num2str(p*(length(var_labels3))) ];
    
    disp([char(var_labels3(v)),' , ', report_string])
    
    fig=figure;
    scatter(var1,var2)
    xlabel([char(var_labels3(v)),',',char(day(1))], 'Interpreter', 'none')
    ylabel([char(var_labels3(v)),',',char(day(2))], 'Interpreter', 'none')
    title(report_string,'Interpreter','none')%Bonferroni-correct for 7 tests (we do not test the interictal FR)

    clearvars r LB UB F df1 df2 p
end

