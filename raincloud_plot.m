function h = raincloud_plot(X,yaxis_val,cl)

% % % [a,b] = ksdensity(X);
% % % 
wdth = 1; % width of boxplot
% % % % TODO, should probably be some percentage of max.height of kernel density plot
% % % 
% % % % density plot
% % % h{1} = area(b,a); hold on
% % % set(h{1}, 'FaceColor', cl);
% % % set(h{1}, 'EdgeColor', [0.1 0.1 0.1]);
% % % set(h{1}, 'LineWidth', 1);
% % % 
% make some space under the density plot for the boxplot
ylim([0 1])
yl = get(gca,'YLim');
set(gca,'YLim',[-yl(2) yl(2)]);

% jitter for raindrops
jit = yaxis_val;


% info for making boxplot
Y = round(quantile(X,[0.25 0.75 0.5]),2);
temp = round(1.5*(Y(2)-Y(1)),2);
n=1;
for i=temp+Y(2):-0.001:Y(2)
    for j=1:size(X,1)
        if round(i,3)==round(X(j),2)
            count(n)=i;
            n=n+1;
            break;
        end
    end
end
Y(4)=count(1);

count=[];
n=1;
for i=Y(1)-temp:0.001:Y(1)
    for j=1:size(X,1)
        if round(i,3)==round(X(j),2)
            count(n)=i;
            n=n+1;
            break;
        end
    end
end
Y(5)=count(1);

% raindrops
h{3} = scatter(X,jit-1,'o','filled');
h{3}.SizeData = 5;
h{3}.MarkerFaceColor = cl;
h{3}.MarkerEdgeColor = cl;

% 'box' of 'boxplot'
h{2} = rectangle('Position',[Y(1) -1-(wdth*yl(2)/2) Y(2)-Y(1) wdth]);
set(h{2},'EdgeColor','k')
set(h{2},'LineWidth',1);
% could also set 'FaceColor' here as Micah does, but I prefer without

% mean line
h{3} = line([Y(3) Y(3)],[-1-(wdth*yl(2)/2) (-1-(wdth*yl(2)/2) +wdth)],'col','k','LineWidth',1);

% whiskers
h{4} = line([Y(2) Y(5)], [(-1-(wdth*yl(2)/2) +wdth/2) (-1-(wdth*yl(2)/2) +wdth/2)],'col','k','LineWidth',1);
h{5} = line([Y(1) Y(4)],[(-1-(wdth*yl(2)/2) +wdth/2) (-1-(wdth*yl(2)/2) +wdth/2)],'col','k','LineWidth',1);
h{6} = line([Y(4) Y(4)],[-1-(wdth*yl(2)/2) (-1-(wdth*yl(2)/2) +wdth)],'col','k','LineWidth',1);
h{7} = line([Y(5) Y(5)],[-1-(wdth*yl(2)/2) (-1-(wdth*yl(2)/2) +wdth)],'col','k','LineWidth',1);

ylim([-4 4])
xlim([0.45 0.9])



