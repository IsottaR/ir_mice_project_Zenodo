% function scatter_box_plot(x,y)
function fig=scatter_box_plot(x,y,boxlabels)

fig=figure
boxplot([x' y'],boxlabels,'Colors',[0 0 0])
for subj=1:length(x)
    
    pos_x=1+(rand(1)-0.5)/10;
    pos_y=2+(rand(1)-0.5)/10;
    hold on
    scatter(pos_x,x(subj),40,'MarkerEdgeColor',[201 198 187]/255,'MarkerFaceColor',[201 198 187]/255,'LineWidth',1.5);
    hold on
    scatter(pos_y,y(subj),40,'MarkerEdgeColor',[201 198 187]/255,'MarkerFaceColor',[201 198 187]/255,'LineWidth',1.5);
    hold on
    
    plot([pos_x pos_y],[x(subj) y(subj)],'k')
    hold on
end
hold off