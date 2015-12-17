function [M,R,C,G]= statistics_1(H,n,m,p,L,D)
% Calculates data relating to several statstics relating to MBL
% M     - Eigendecomposition of the matrix
% R     - Level spacing ratio 
% C     - Infinite time 2-site pauli correlators
% G     - Thoulless conductivity
% recquires input of
% H     - Hamiltonian 
% n,m,p - local basis on operators
% L     - length of system
% D     - L-by-L distance matrix

% system dim
dim=2^L;

% check dimensions
assert(size(H,1)==dim);
assert(size(H,2)==dim);
for i=1:L
    assert(size(n{i},1)==dim);
    assert(size(n{i},2)==dim);
    assert(size(p{i},1)==dim);
    assert(size(p{i},2)==dim);
    assert(size(m{i},1)==dim);
    assert(size(m{i},2)==dim);
end



% identity
I=speye(dim);


% find the eigenvalues of H = U E U'
[U,E]=eig(full(H));
% sort them
[~,J]=sort(diag(E));
U=U(:,J);
E=E(J,J);
% get diagonal elements
E_d=diag(E);

M.U=U;
M.E=E_d;

% PART 1 - LEVEL SPACING

% find level diffs
dE=E_d(2:end)-E_d(1:end-1);
% find ratios
rk=dE(1:end-1)./dE(2:end);
% take min of rk and 1/rk so that rk in [0,1]
rk=min(rk,rk.^(-1));
% take average level spacing ratio
mr=mean(rk);

R.rk=rk;
R.mr=mr;

% PART 2 - Infinite time 2-site pauli correlators

% initialise correlations
Ci=zeros(L,L);

% initialise pauli diagonal element matrices
S0=cell(1,3);
S0(:)={zeros(dim,L)};
S1=cell(1,3);
S1(:)={zeros(dim-1,L)};
% for each site i evaluate <n|\sigma_i|n> for each pauli \sigma and each
% eigenvector |n>
for i=1:L
    % create the three pauli matrices
    s{1}=p{i}+m{i};
    s{2}=1i*(m{i}-p{i});
    s{3}=2*n{i}-I;
    % convert to the eigenbasis and get diagonal elements
    for k=1:3
        % get diagonal in the eigenbasis
        for a=1:dim
            V=s{k}*U(:,a);
            S0{k}(a,i)=(U(:,a)')*V;
            if a==dim
                break
            end
            S1{k}(a,i)=(U(:,a+1)')*V;
        end
    end
end
% for each pair of sites find the maximal correlation
for i=1:L
    for j=(i):L
        % make a 3x3 matrix of the correlators
        cij=zeros(3,3);
        % fill in the elements
        for k=1:3
            for m=1:3
                cij(k,m)=dot(S0{k}(:,i),S0{m}(:,j))/dim;
            end
        end
        % find the max singular value of cij
        Ci(i,j)=max(svd(cij));
        Ci(j,i)=Ci(i,j);
    end
end
% number of pairs
np=L*(L+1)/2;
% store log correlations between pairs with distance
V=zeros(np,2);
k=1;
for i=1:L
    for j=i:L
        V(k,1)=D(i,j);
        V(k,2)=log(Ci(i,j));
        k=k+1;
    end
end
% theil-sen line fit gives the inverse localisation length
[b,db]=theil_sen_slope_fit(V(:,1),-V(:,2));

C.C=Ci;
C.V=V;
C.b=b;
C.db=db;

% PART 3 thouless conductances

% initialise data
Gk=cell(1,3);
Gk(:)={zeros(dim-1,L)};

% for each site try with V= sigma_i for each i
for k=1:3
    for i=1:L
        % calculate G=log(V_{n+1,n}/(E_{n+1}-E_n))
        Gk{k}(:,i)=log(abs(S1{k}(:,i))./(dE));
    end
end
% get the mean values for the three perturbatiosn
mG=zeros(1,3);
for k=1:3
    mG(k)=median(Gk{k}(:));
end
% energy level means
mE=(E_d(2:end)+E_d(1:end-1))/2;

G.Gk=Gk;
G.mG=mG;
G.mE=mE;

end

