function emitters_filt = filterEmitters_end(emitters)

keep = false(1, numel(emitters));

for k = 1:numel(emitters)
    if ~emitters(k).alive
        keep(k) = true;
    end
end

emitters_filt = emitters(keep);
num = numel(emitters);
num_filt = numel(emitters_filt);
filtered = (num-num_filt)/num*100;
fprintf('Totally %d emitters left. %.2f%% end-frame-alive emitters are filtered.\n', num_filt, filtered);
end