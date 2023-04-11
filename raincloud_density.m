function h = raincloud_density(X,cl)

[a,b] = ksdensity(X);

wdth = 0.8; % width of boxplot
% TODO, should probably be some percentage of max.height of kernel density plot

% density plot
h{1} = area(b,a); hold on
set(h{1}, 'FaceColor', cl);
set(h{1}, 'EdgeColor', [0.1 0.1 0.1]);
set(h{1}, 'LineWidth', 1);