function []=zip_all_files(zip_fname,fname_pattern)
% zips all the files with a certain filename pattern in their current
% directory

[dir,~,~]=fileparts(fname_pattern);

zip(fullfile(dir,zip_fname),all_file_names(fname_pattern));

end

