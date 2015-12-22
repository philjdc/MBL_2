function [] = data_gen_east_model_stats(params,run_name,pool_size)

% get date_time
date_time=strrep(strcat(datestr(clock,'yyyy-mm-dd-HH--MM'),'m',datestr(clock,'ss'),'s'),'--','h');
% number of instances
N=length(params);
% make a vector to label which instances are complete
complete=zeros(1,N);
% initialise all tasks as incomplete
complete(:)=false;  

% initialise complete filename
complete_fname='complete_tmp';

% if output dir doesnt exist make it
out_dir=fullfile(pwd,'out');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

% if tmp dir doesn't exist make it then and go there, if already
% existing look for tmp files from existing run
tmp_dir=fullfile(pwd,'tmp',run_name);

if ~exist(tmp_dir,'dir')
    mkdir(tmp_dir);
    wor_dir=cd(tmp_dir);
    % make a complete_tmp file 
    save(complete_fname,'complete');
else
    wor_dir=cd(tmp_dir);
    if exist(strcat(complete_fname,'.mat'),'file')
        load(strcat(complete_fname,'.mat'));
        fprintf(strcat('Found existing partial execution of "',run_name,'". Continuing execution.\n'));
    end
end

% open matlab pool
if pool_size>1
    matlabpool('open',pool_size);
end

% define an output pattern
data_fname_pattern=strcat(run_name,'_*');

% start collecting data
parfor n=1:N
    % if not complete
    if ~complete(n)
        % make a filename
        data_fname=strrep(data_fname_pattern,'*',num2str(n));
        % assert that it doesnt exist
        assert(~exist(data_fname,'file'));
        % get params
        p=params{n};
        % collect data
        data=statistics_1_handler_east_model(p.J,p.t,p.h,p.dh,p.L,p.left_flip,p.right_flip,p.boundary_cyclic,p.N);
        %data=p;
        % make filename and save data
        save_data(data_fname,data);
    end
    % set as complete
    complete(n)=true;
    % record completeion
    save_complete(complete_fname,'complete');
    % print progress
    
    fprintf(['Case ',num2str(n),' of ',num2str(N),' complete.\n']);
end

% close pool
if pool_size>1
    matlabpool('close');
end

% zip all the files
zip_fname=strcat(run_name,'_',date_time);
zip_all_files(zip_fname,fullfile(pwd,data_fname_pattern));
zip_fname=fullfile(tmp_dir,strcat(zip_fname,'.zip'));

% go back to working dir
cd(wor_dir);

% move the tmp_dir to out
fclose('all');
movefile(zip_fname,out_dir);

% delete temporary files
rmdir(tmp_dir,'s');

end

function []=save_complete(fname,complete)
% hides save function from par pool

save(fname,'complete');

end

function []=save_data(fname,data)
% hides save function from par pool

save(fname,'data');

end

