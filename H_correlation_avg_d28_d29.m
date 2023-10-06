%Correlate AVG(GE28,GE29) with AVG(number of seizure)

% - AVERAGE number of seizures
% - correlate average graph measures with average number of seizures

% -------> change path at line 19

clc
clear
close all

%% subject initialisation
ses_id='ses-d28-d29'; %time points of observation
subj=[12 13 15:20 22:33]; %animal ID
%--------------------------------------

%% variable initialisation

BIDSfolder='H:\Isotta\DATA\ir_mice_project\RS\data2publish';
task='task-rest';%
f_band_name=[ {'delta'}, {'lowTheta'},{'highTheta'} ,{'beta'} ,{'gamma'}, {'broadband'}];
band2analyse=[2 3];
var_labels={'LI','GE_R','CC_R'};%only the metrics (glob and hemi) that were significant

freq_band_freq=[1 4; 4 8; 8 12; 12 30; 30 40; 1 40];

%% load EA info
EAinfo28=readtable(fullfile(BIDSfolder,'EA_info.xlsx'),'Sheet','d28');
EAinfo29=readtable(fullfile(BIDSfolder,'EA_info.xlsx'),'Sheet','d29');
var2_labels=EAinfo28.Properties.VariableNames(2:8);

%% load network measures at day 28
for s=1:length(subj)
    
    clearvars GE SO SO_h LI PDC
    
    %subject id
    sub_id=subj(s);
    
    %---------------------- Load network data -------------------------
    load(fullfile(BIDSfolder,'derivatives','network_metrics_wpli',['sub-',sprintf('%02d',sub_id)],ses_id,'eeg',...
        ['sub-',sprintf('%02d',sub_id),'_',ses_id,'_',task,'_network_metrics.mat']));
    
    %-----------global measures
    GE_2corr(s,:)=GE;
    LI_2corr(s,:)=LI;
    avgCC_2corr(s,:)=avgCC;
    
    %-----------hemispheric measures
    
    GE_R_2corr(s,:)=GE_hemi(1,:);
    GE_L_2corr(s,:)=GE_hemi(2,:);
    
    CC_R_2corr(s,:)=CC_hemi(1,:);
    CC_L_2corr(s,:)=CC_hemi(2,:);


    %------------ get behavioral measures of each animal --------------
    sub_idx28=find(table2array(EAinfo28(:,1))==sub_id);
    measures28(s,:)=table2array(EAinfo28(sub_idx28,2:8));
    
    sub_idx29=find(table2array(EAinfo29(:,1))==sub_id);
    measures29(s,:)=table2array(EAinfo29(sub_idx29,2:8));
    
    subjects28(s,1)=table2array(EAinfo28(sub_idx28,1));
    subjects29(s,1)=table2array(EAinfo29(sub_idx29,1));
    
    if (subjects28(s,1)~=subj(1,s)) || (subjects29(s,1)~=subj(1,s))
        error('the subjects are not consistent')
    end
end

%% correlation
%correlation between global efficiency (in all freq band) and numb of
%seizures
for b=1:length(band2analyse)
    for v=1:length(var_labels)
        for m=1:size(var2_labels,2)
            
            %network metric
            eval(['var1=squeeze(',char(var_labels(v)),'_2corr(:,band2analyse(b)));']);
            
             %---------------------- (metric28+metric29)/2 ----------------------
            %epilepsy metrics
            var2_28=measures28(:,m);
            var2_29=measures29(:,m);
            var2=(var2_28+var2_29)./2;
            
            if any(isoutlier(var1))||any(isoutlier(var2))
                corrtype='Spearman';
            else
                corrtype='Pearson';
            end
            
            [rho,pval]=corr(var1,var2,'Type',corrtype)

            if pval<.05
  
                fig=figure('Name',[corrtype,' corr cofficient in ',char(f_band_name(band2analyse(b)))])
                scatter(var1,var2)
                hold on
                fitL=polyfit(var1,var2,1);
                plot(var1,polyval(fitL,var1),'k')
                hold off
                xlabel([char(var_labels(v)),'  in ',char(f_band_name(band2analyse(b)))])

                ylabel(char(var2_labels{1, m}))
                title(['rho: ',num2str(rho),' p:',num2str(pval)])
            end
        end
    end
end


%% scatter plot (according to groups of animals) --> Fig 6

% --------------- find which animals had seizures and define groups -------
m=4 % plot it for seizures only
%epilepsy metrics
var2_28=measures28(:,m);
var2_29=measures29(:,m);

%find mice who had zero seizures
seizures_28=find(var2_28==0);
seizures_29=find(var2_29==0);

%mice who had no seizures on neither days
seizure28_29=intersect(seizures_28,seizures_29);

grp=ones(size(var2_28))*3; %seizures on both days=3
grp(seizures_28)=2;% no seizures on day 28 =2
grp(seizure28_29)=1;% no seizures on neither days =1;

% --------------- organize and plot the data (LI) -------------------------
var_labels={'LI'}
for b=1:length(band2analyse)
 
    for v=1:length(var_labels)
        %network metric
        eval(['var1(:,b)=squeeze(',char(var_labels(v)),'_2corr(:,band2analyse(b)));']);
        lab{1,b}=[char(var_labels(v)),'_',char(f_band_name(band2analyse(b)))];
    end
end
    
fig=figure
boxplot(var1)
for v=1:size(var1,2)
    hold on
    gscatter(ones(length(grp),1).*(v+(rand(length(grp),1)-0.5)/10),var1(:,v),grp,[[187 204 51]/255 ;[238 221 136]/255;[238 136 102]/255],'.',20)
    legend('No iHPD','iHPD on one day','iHPD on both days')
end
xticklabels(lab)

% --------------- organize and plot the data (HE and HCC) -----------------
%organize
clear var_labels
b=2
var_labels={'GE_R','CC_R'}

for v=1:length(var_labels)
    %network metric
    eval(['var1(:,v)=squeeze(',char(var_labels(v)),'_2corr(:,band2analyse(b)));']);
    lab{1,v}=[char(var_labels(v)),'_',char(f_band_name(band2analyse(b)))];
end

%actually plot  
fig=figure
boxplot(var1)
for v=1:size(var1,2)
    hold on
    gscatter(ones(length(grp),1).*(v+(rand(length(grp),1)-0.5)/10),var1(:,v),grp,[[187 204 51]/255 ;[238 221 136]/255;[238 136 102]/255],'.',20)
    legend('No iHPD','iHPD on one day','iHPD on both days')
end
xticklabels(lab)
