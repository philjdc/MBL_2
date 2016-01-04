function agg_data=data_aggregate(varargin)
% aggregates data from a run 


% set up directories
% data source (output directory)
out_dir=fullfile(pwd,'out');
% working directory
tmp_dir=fullfile(pwd,'tmp');
% unzip to a temp forlder in working dir
uzp_dir=fullfile(tmp_dir,'agg');
mkdir(uzp_dir);

% get the file name
if nargin==0
    [fname,fdir]=uigetfile(fullfile(out_dir,'*.zip'));
    fname=fullfile(fdir,fname);
else
    fname=varargin{1};
end

% unzip the files
file_list=unzip(fname,uzp_dir);
% count them
N=length(file_list);

% initialise aggregate data
agg_data.dh=zeros(1,N);
agg_data.J=zeros(1,N);
agg_data.t=zeros(1,N);
agg_data.h=zeros(1,N);
agg_data.L=zeros(1,N);
agg_data.r=zeros(1,N);
agg_data.rm=zeros(1,N);
agg_data.b=zeros(1,N);
agg_data.db=zeros(1,N);
agg_data.gm1=zeros(1,N);
agg_data.gm2=zeros(1,N);
agg_data.gm3=zeros(1,N);
agg_data.g1=zeros(50,N);
agg_data.g2=zeros(50,N);
agg_data.g3=zeros(50,N);


% iterate through them
for n=1:N
    % get the current file
    data=load(file_list{n});
    data=data.data;
    
    agg_data.dh(n)=data.P.dh;
    agg_data.J(n)=data.P.J;
    agg_data.t(n)=data.P.t;
    agg_data.h(n)=data.P.h;
    agg_data.L(n)=data.P.L;
    agg_data.r(n)=data.R.C.mean;
    agg_data.rm(n)=data.R.C.med;
    agg_data.b(n)=data.C.b;
    agg_data.db(n)=data.C.db;
    agg_data.gm1(n)=data.G{1}.C.mean;
    agg_data.gm2(n)=data.G{2}.C.mean;
    agg_data.gm3(n)=data.G{3}.C.mean;
    agg_data.g1(:,n)=data.G{1}.W1.mE;
    agg_data.g2(:,n)=data.G{2}.W1.mE;
    agg_data.g3(:,n)=data.G{3}.W1.mE;
    
end

% delete the tmp files
rmdir(uzp_dir,'s');


end

