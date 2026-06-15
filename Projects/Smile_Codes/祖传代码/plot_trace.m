function plot_trace(f_series, ex_time, out_put_dir)

[loc_num, frames] = size(f_series);
time = linspace(0, frames-1, frames) * ex_time;

figure;
set (gcf,'Position',[300 100 600 500], 'color','w');
set(0,'DefaultFigureVisible', 'off');
for num = 1:loc_num
    I = plot(time, f_series(num,:), 'black', 'LineWidth', 2);
    title(['Si-Rhodamine3: trace', num2str(num), ' (0.29 kW cm^{-2})']);
    xlabel('Time (s)');
    ylabel('Photons per frame');
    xlim([0,100]);
    set(gca, 'FontName', 'Times New Roman');
    set(findall(gcf,'-property','FontSize'), 'FontSize', 24);
    
    saveas(I, [out_put_dir, 'trace', int2str(num), '.tif']);
end

end

