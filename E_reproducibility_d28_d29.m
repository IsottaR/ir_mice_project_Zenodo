%Investigate the level of variability in adjecency matrices (all frequency
%bands) between day 28 and day 29
clear all
close all
clc

%% subject initialisation
day={'d28' 'd29'}; %time points of observation
subj=[12 13 15:20 22:33]; %animal ID
%--------------------------------------

%% variable initialisation
BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\data2publish\derivatives';
task='task-rest';
derivative_folder='wpli';
f_band_name=[ {'delta'}, {'lowTheta'},{'highTheta'} ,{'beta'} ,{'gamma'}, {'broadband'}];

%% load EA info
EAinfo28=readtable('H:\Isotta\DATA\ir_mice_project\RS\data2publish\EA_info.xlsx','Sheet','d28');
EAinfo29=readtable('H:\Isotta\DATA\ir_mice_project\RS\data2publish\EA_info.xlsx','Sheet','d29');
var_labels3=EAinfo28.Properties.VariableNames(2:8);

%% organize epilepsy parameters
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

%% final folder
final_folder=fullfile(BIDSfolder,['stats_',derivative_folder,'_reproducibility_d28d29']);
if ~exist(final_folder)
   mkdir(final_folder) 
end

%% organize connectomes and store them
for d=1:length(day)
    %get network metric for each session: d0, d28 or d29
    %session ID
    ses_id=['ses-',char(day(d))];
    for s=1:length(subj)
        clearvars GE STR LI CC
        
        %subject id
        sub_id=subj(s);
        
        %load connectivity matrix
        load(fullfile(BIDSfolder,derivative_folder,['sub-',sprintf('%02d',sub_id)],ses_id,'eeg',...
            ['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,'.mat']));
        
        %save metrics
        wpli_all(:,:,:,d,s)=cat(3,wpli_avg_delta, wpli_avg_lowTheta,...
            wpli_avg_highTheta,wpli_avg_beta,wpli_avg_gamma,wpli_avg);
    end
end


%% do correlation analyses on CONNECTION VALUES (functional connectomes) -> make Fig 5.a
for b=1:size(wpli_all,3)
    for s=1:size(wpli_all,5)
        clear curr_mat_x curr_mat_y var 1 var2
        curr_mat_x=squeeze(wpli_all(:,:,b,1,s));%day28
        curr_mat_y=squeeze(wpli_all(:,:,b,2,s));%day29
        
        %vectorize upper matrix
        var1=curr_mat_x(triu(true(size(curr_mat_x)),1));
        var2=curr_mat_y(triu(true(size(curr_mat_y)),1));
        
        [rho(s,b),pval(s,b)]=corr(var1,var2,'Type','Pearson')
    end
    
    %average across animals
    clear curr_mat_x curr_mat_y var 1 var2

    curr_mat_x=mean(squeeze(wpli_all(:,:,b,1,:)),3);%day28
    curr_mat_y=mean(squeeze(wpli_all(:,:,b,2,:)),3);%day29
    
    %vectorize upper matrix
    var1=curr_mat_x(triu(true(size(curr_mat_x)),1))
    var2=curr_mat_y(triu(true(size(curr_mat_y)),1))
    
    [rho_avg(b),pval_avg(b)]=corr(var1,var2,'Type','Pearson')
    
    fig=figure
    scatter(var1,var2)
    hold on
    fitL=polyfit(var1,var2,1);
    plot(var1,polyval(fitL,var1),'k')
    hold off
    xlabel(['avg ',char(f_band_name(b)),' FC, ',char(day(1))])
    
    ylabel(['avg ',char(f_band_name(b)),' FC, ',char(day(2))])
    title(['rho: ',num2str(rho_avg(b)),' p:',num2str(pval_avg(b)*5)],'Interpreter','none') %Bonferroni-correct
end

%% do correlation analyses on epilepsy PARAMETERS (saeizures etc.) -> make rest of Fig 5

for v=1:length(var_labels3)
    %---------------------- (metric28+metric29)/2 ----------------------
    eval(['var1=squeeze(measures2plot(1,:,v));']);%day 28
    eval(['var2=squeeze(measures2plot(2,:,v));']);%day 29
    
    [rho_avg,pval_avg]=corr(var1',var2','Type','Pearson');
    
    %try calculating the intra-class correlation coefficient
    [r, LB, UB, F, df1, df2, p] = ICC([var1; var2]', 'C-1');
    
    fig=figure
    scatter(var1,var2)
    hold on
    fitL=polyfit(var1,var2,1);
    plot(var1,polyval(fitL,var1),'k')
    hold off
    xlabel(char(day(1)), 'Interpreter', 'none')
    ylabel(char(day(2)), 'Interpreter', 'none')
    title([char(var_labels3(v)), ' rho: ',num2str(rho_avg),' p:',num2str(pval_avg*(length(var_labels3)-1))],'Interpreter','none')%Bonferroni-correct for 7 tests (we do not test the interictal FR)
end

