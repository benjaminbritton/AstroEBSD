function [pattern_list,num_patterns]=Folder_Prep(InputUser)
    


if exist(InputUser.EBSD_folder,'dir') ~=7
error(['The folder ' InputUser.EBSD_folder ' does not exist']);
end

folder_dir = dir(InputUser.EBSD_folder);
[~,idx] = sort([folder_dir.datenum]);
folder_dir=folder_dir(idx);
folder_files={folder_dir.name};
%find the bmp names
bmp_names=folder_files(logical(1-cellfun('isempty',strfind(folder_files,'bmp'))));
%find the tif names
tif_names=folder_files(logical(1-cellfun('isempty',strfind(folder_files,'tif'))));

all_patterns=[bmp_names tif_names];
num_patterns=numel(all_patterns);

pattern_list=cell(num_patterns,1);
for n=1:num_patterns
    pattern_list{n}=fullfile(InputUser.EBSD_folder,all_patterns{n});
end

end