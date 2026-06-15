function pos_ff = checkFrameEmitters(pos_matrix, data, ONframes)
% Visualize emitters active within chosen frames.
% [ONframes] may be a list.

    % --- mask of emitters active within first ONframe ---
    on_mask = any(~isnan(pos_matrix(ONframes,1,:)), 1);
    on_mask = squeeze(on_mask);

    % --- mean position ---
    pos_tmp = pos_matrix(ONframes,:,on_mask);
    pos_ff  = squeeze(mean(pos_tmp,1,'omitnan') + 0.5).';
    pos_ff  = round(pos_ff);

    % --- visualization ---
    figure
    viz.plotImage(data,1:10,'gray', ...
        sprintf('Emitters active within the chosen %d frames', numel(ONframes)));

    hold on
    scatter(pos_ff(:,2), pos_ff(:,1), 40, 'g', 'x')

    fprintf('Total %d emitters active within the chosen %d frames.\n', ...
            sum(on_mask), numel(ONframes));
end