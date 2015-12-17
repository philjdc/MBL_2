function fnames = all_file_names(fname_pattern)
% given a string fname_pattern which is a full pathname containg the 
% possibly wild card charachter *, this return all files that match the
% pattern.

% get file name parts
[new_dir,name,ext]=fileparts(fname_pattern);
fname=strcat(name,ext);
% move to directoy
old_dir=cd(new_dir);
% get list of matching filenames
contents=dir(fname);
N=length(contents);
fnames=cell(N,1);
for n=1:N
    fnames{n}=contents(n).name;
end
% move back
cd(old_dir);


end

