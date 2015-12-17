function []=include_all()
% adds all directorites in current working directory to the path

% get current dir
curr_dir=pwd;
addpath(curr_dir);
fprintf('addpath: %s \n',curr_dir);
% get list of folder contents
conts=dir;
for i=1:length(conts)
    if (~(conts(i).name(1)=='.'))
        % make filename
        path_i=fullfile(curr_dir,conts(i).name);
        if isdir(path_i)
            addpath(path_i);
            fprintf('addpath: %s \n',path_i);
            cd(path_i);
            include_all;
            cd(curr_dir);
        end
    end
end
