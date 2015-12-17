function Data = statistics_1_handler_east_model(J,t,h,dh,L,left_flip,right_flip,boundary_cyclic,N)
% amalgamates statistics for N copies of the east model hamiltonian
%
% Arguments are the parameters of the East Model hamiltonian, with N which
% is the number of Hamiltonians to sample.
%
% Output is object Data with fields
% Data
%     .P    -   parameters of the data collection 
% 	  .R    -   Level statistics
%       .rk -   all level spacing ratios
%       .mr -   mean level spacing ratio
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
R.rk=[];
C.C=cell(1,N);
C.V=[];
G.Gk=cell(1,3);
G.Gk(:)={[]};
G.mE=[];
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
            flag_success=false;
        end
        if flag_success
            % if successfull continue
            break;
        end
    end
    
    % amalgamate data
    E=[E;Mn.E];
    R.rk=[R.rk;Rn.rk];
    C.C{nn}=Cn.C;
    C.V=[C.V;Cn.V];
    for k=1:3
        G.Gk{k}=[G.Gk{k};Gn.Gk{k}];
    end
    G.mE=[G.mE;Gn.mE];
    
    
    
    
    
    
end

% calculate centres and widths of some of the distributions to aid analysis
% mean level spacing ratio
R.rm=mean(R.rk);
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
% estimate inverse loc. length
[b,db]=theil_sen_slope_fit(C.V(:,1),-C.V(:,2));
C.b=b;
C.db=db;
% mean thoulless conductances
G.mG=zeros(1,3);
G.dG=zeros(1,3);
for k=1:3
    G.mG(k)=median(G.Gk{k}(:));
    G.dG(k)=median(abs(G.Gk{k}(:)-G.mG(k)));
end
% rescale energies to [0,1]
E_max=max(E);
E_min=min(E);
G.mE=(G.mE-E_min)/(E_max-E_min);

% output data
Data.P=P;
Data.R=R;
Data.C=C;
Data.G=G;

end

