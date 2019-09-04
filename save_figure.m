function save_figure(figure_name,fig_path,save_fig,save_png,save_eps)
%
% Syntax: save_figure(figure_name,fig_path,save_fig,save_png,save_eps)
%
% This function gets a path a figure name and saves is in any of the fig,
% png and eps formats in that path. 
%
% Input: save_fig, save_png and save_eps are flags that determine which format
%        should be save. 
%
% Created by Mehrdad Pourfathi on 3/3/2014

if save_png
    print(gcf,'-dpng','-r300',[fig_path figure_name '.png'])
end
if save_eps
    saveas(gcf,[fig_path figure_name '.eps'],'eps2c');
end
if save_fig
    saveas(gcf,[fig_path figure_name '.fig'],'fig');
end 