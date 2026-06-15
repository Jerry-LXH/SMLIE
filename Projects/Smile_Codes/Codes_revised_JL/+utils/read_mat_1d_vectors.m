function data_all = read_mat_1d_vectors(dataDir, pattern)

files = dir(fullfile(dataDir, pattern));
data_all = [];

for i = 1:numel(files)
    if files(i).isdir || startsWith(files(i).name, '.') || startsWith(files(i).name, '._')
        continue;
    end

    filePath = fullfile(files(i).folder, files(i).name);
    tmp = load(filePath);
    fn = fieldnames(tmp);
    x = tmp.(fn{1});
    data_all = [data_all; x(:)];
end

end