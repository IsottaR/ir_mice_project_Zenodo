function bad_ch=find_bad_channels_mice(data)

%remove channel 31 (should be zero)
cfg=[];
cfg.channel={'all', '-31'};
data=ft_selectdata(cfg,data);

%get SD of channels over epochs
std_ch=std(cat(2,data.trial{:}), [], 2);

bln = isoutlier(sqrt(std_ch), 'quartiles', 'ThresholdFactor',2); %HightTh=Q3+2*(Q3-Q1)
% remove ref from bad channels (std of Ref = 0)
bln = bln & std_ch>0;
bad_ch = data.label(bln);

if isempty(bad_ch)
    bad_ch = {};
end

end
