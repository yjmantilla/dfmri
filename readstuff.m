clear all
data_path = 'Z:\dfmri\SAN\Rest';
MyFolderInfo = dir(data_path);
files = {MyFolderInfo.name}; % i hate matlab... 
file_names = {};
j = 1;
bytes = {MyFolderInfo.bytes};
for i = 1:length(files)
    if bytes{i}> 0
        file_names{j} = files{i};
        j = j + 1;
    end
end