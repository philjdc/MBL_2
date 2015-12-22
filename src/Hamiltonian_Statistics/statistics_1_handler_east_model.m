function Data = statistics_1_handler_east_model(J,t,h,dh,L,left_flip,right_flip,boundary_cyclic,N)
% amalgamates statistics for N copies of the east model hamiltonian
%
% Arguments are the parameters of the East Model hamiltonian, with N which
% is the number of Hamiltonians to sample.
%
% Output is object Data with fields
% Data
%     .P    -   parameters of the data collection 
% 	  .R    -   Level statistics distribution data
%     .C    -   Correlation statistics
%       .C  -   correlations between single site operators averaged over
%               time, and maximised over operators
%       .V  -   V(:,1) distance between sites, V(:,2) log correlation
%               between sites
%       .b  -   inverse localisation length obtained from theil-sen linear
%               fit of the data in V
%       .db -   error on db (median absolute deviation of all estimates of
%               b from the median estimate.
%     .G        Thouless conductivity statistics
%       .Gk(i)  Thouless conductances obtained from local perturbations of
%               sigma(i) (pauli matrices)
%       .mG(i)  median thouless conductance for each of the Gk(i)
%       .dG(i)  median absolute deviation of the Gk(i) from mG(i)
%       .mE     the average of the two energy levels used to calculate each
%               Thouless conductance, allowing for anaylsis of possible
%               mobility edges

% standard vals
nWindow=50;
nCoeffs=30;


% get date_time
date_time=strrep(strcat(datestr(clock,'yyyy-mm-dd-HH--MM'),'m',datestr(clock,'ss'),'s'),'--','h');

% store params
P.J=J;
P.t=t;
P.h=h;
P.dh=dh;
P.L=L;
P.left_flip=left_flip;
P.right_flip=right_flip;
P.boundary_cyclic=boundary_cyclic;
P.N=N;
P.date_time=date_time;

% initialise data
rk=[];
C.C=cell(1,N);
V=[];
Gk=cell(1,3);
Gk(:)={[]};
mE=[];
E=[];

% get the zero field hamiltonian
[H0,n,p,m,~,~,D]=hamiltonian_east_model(J,t,h,0,L,left_flip,right_flip,boundary_cyclic);

% iterate
for nn=1:N


    
    % try disorder the hamiltonian and get the statistics
    while true
        % flag for retry
        flag_success=true;
        % initialise the hamiltonian
        H=H0;
        try
            % disorder the hamiltonian
            for i=1:L
                H=H+n{i}*unifrnd(-dh/2,dh/2);
            end
            % get statistics
            
            [Mn,Rn,Cn,Gn]= statistics_1(H,n,m,p,L,D);
        catch
            % catch error and flag it
        %    flag_success=false;
        end
        if flag_success
            % if successfull continue
            break;
        end
    end
    
    % amalgamate data
    E=[E;Mn.E];
    rk=[rk;Rn.rk];
    C.C{nn}=Cn.C;
    V=[V;Cn.V];
    for k=1:3
        Gk{k}=[Gk{k};Gn.Gk{k}];
    end
    mE=[mE;Gn.mE];
    
    
    
    
    
    
end

% calculate centres and widths of some of the distributions to aid analysis
% mean level spacing ratio
R=sample_data(rk,nCoeffs,[0,1]);
% median correlations
C.mC=zeros(L,L,N);
for nn=1:N
    C.mC(:,:,nn)=C.C{nn};
end
C.dC=C.mC;
C.mC=median(C.mC,3);
for nn=1:N
    C.dC(:,:,nn)=abs(C.dC(:,:,nn)-C.mC);
end
C.dC=median(C.dC,3);
% collate data too
dVals=unique(V(:,1));
C.mV=zeros(1,max(dVals)+1);
C.dV=zeros(1,max(dVals)+1);
for d=1:max(dVals)
    % get all cvals at this distance
    I=(V(:,1)==d);
    C.mV(d+1)=median(V(I,2));
    C.dV(d+1)=median(abs(V(I,2)-C.mV(d+1)));
end
% estimate inverse loc. length
[b,db]=theil_sen_slope_fit(V(:,1),-V(:,2));
C.b=b;
C.db=db;
% mean thoulless conductances

% rescale energies to [0,1]
E_max=max(E);
E_min=min(E);
mE=(mE-E_min)/(E_max-E_min);
% find median G in each window
G=cell(1,3);

W1mE=zeros(nWindow,1);
W1dE=zeros(nWindow,1);
W1BE=zeros(nWindow,1);
W1N=zeros(nWindow,1);
W2mE=zeros(nWindow,1);
W2dE=zeros(nWindow,1);
W2N=zeros(nWindow,1);
W2BE=zeros(nWindow,1);
for k=1:3
    G{k}=sample_data(Gk{k}(:),nCoeffs);
    W1BE=(0:nWindow)'/nWindow;
    W2BE=[0,quantile(mE,nWindow)]';
    for i=1:nWindow
        % get evenly spaced windows
        I=(mE>W1BE(i))&(mE<=W1BE(i+1));
        G_window=Gk{k}(I,:);
        G_window=G_window(:);
        % get data for evenly spaced windows
        W1mE(i)=median(G_window);
        W1dE(i)=median(abs(G_window-W1mE(i)));
        W1N(i)=length(G_window(:));
        
        % get evenly populated windows
        I=(mE>W2BE(i))&(mE<=W2BE(i+1));
        G_window=Gk{k}(I,:);
        G_window=G_window(:);
        % get data for evenly populated windows
        W2mE(i)=median(G_window);
        W2dE(i)=median(abs(G_window-W2mE(i)));
        W2N(i)=length(G_window);
    end
    G{k}.W1.mE=W1mE;
    G{k}.W1.dE=W1dE;
    G{k}.W1.N=W1N;
    G{k}.W1.edges=W1BE;
    G{k}.W2.mE=W2mE;
    G{k}.W2.dE=W2dE;
    G{k}.W2.N=W2N;
    G{k}.W2.edges=W2BE;
end



% output data
Data.P=P;
Data.R=R;
Data.C=C;
Data.G=G;

end

