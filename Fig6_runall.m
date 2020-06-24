% Make sure Figure 1 is called and the docked window is maximized
% set(0,'DefaultFigureWindowStyle','docked');figure(1);
% Generate all the subplots seen in Figure 6
clear;clc;
diary off;
diary('figs/Figure6_console.txt');
diary on;
fprintf('\nNew experiment\n\n%s\n', datestr(now, 'DD mmmm, YYYY HH:MM:SS'));

% Make sure Figure window is open and maximized
for figure_number=[1,2,3]
    Fig6_alg1_vs_alg2;
    set(gcf, 'InvertHardCopy', 'off');
    saveas(gcf, sprintf('figs/Figure6%c.png',figure_number+64), 'png');
    savefig(gcf, sprintf('figs/Figure6%c.fig',figure_number+64), 'compact');
end
diary off;