function check_localization(data_check, loc_check)

figure('visible','on')
imagesc(data_check);
set(gca,'YDir','normal')
% axis square
% title('sample');
% colormap hot  
colorbar
% set(gca, 'FontName', 'Arial');
% set(findall(gcf,'-property','FontSize'),'FontSize',28);
% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
hold on
plot(loc_check(:,2),loc_check(:,1),'o','MarkerSize',9,'Color',[0.8500 0.3250 0.0980]);
set(gcf,'position',[300 100 600 500]);
end

