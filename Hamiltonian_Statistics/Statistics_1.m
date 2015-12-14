function Data = Statistics_1(H,n,m,p,L,D,raw)
% Calculates several statstics relating to MBL
% - Level spacing ratio 
% - Infinite time 2-site pauli correlators
% - Thoulless conductivity

% check dimensions
assert(size(H,1)==2^L);
assert(size(H,2)==2^L);
assert(size(n)==2^L);
assert(size(m)==2^L);
assert(size(p)==2^L);

% multiplicative identity
I=speye(2^L);
% additive identity
Z=sparse(2^L,2^L);

% find the eigenvalues of H = U E U'
[U,E]=eig(full(H));

% store raw data
if raw
    Data.raw.U=U;
    Data.raw.E=E;
    Data.raw.H=H;
end

% PART 1 - LEVEL SPACING

% get eigenvalues
E_d=diag(E);
% find level diffs
dE=E_d(2:end)-E_d(1:end-1);
% find ratios
rk=dE(1:end-1)./dE(2:end);
% take min of rk and 1/rk so that rk in [0,1]
rk=min(rk,rk.^(-1));
% take average level spacing ratio
r1=mean(rk);
% take median level spacing ration
r2=median(rk);

Data.stats.r=r1;

% PART 2 - Infinite time 2-site pauli correlators

% initialise correlations
C=zeros(L,L);

% initialise pauli diagonal element matrices
S=cell(1,3);
for n=1:3
    S{n}=zeros(2^L,L);
end

% for each site i evaluate <n|\sigma_i|n> for each pauli \sigma and each
% eigenvector |n>
for i=1:L
    % create the three pauli matrices
    s{1}=m{i}+p{i};
    s{2}=1i*(p{i}-m{i});
    s{3}=2*m{i}-I;
    % convert to the eigenbasis and get diagonal elements
    for n=1:3
        S{n}(:,i)=diag(U'*s{n}*U);
    end
end
% for each pair of sites find the maximal correlation
for i=1:L
    for j=(i+1):L
        % make a 3x3 matrix of the correlators
        cij=zeros(3,3);
        % fill in the elements
        for n=1:3
            for m=1:3
                cij{n,m}=dot(S{n}(:,i),S{m}(:,j));
            end
        end
        % find the max singular value of cij
        C(i,j)=max(svd(cij));
        C(j,i)=C(i,j);
    end
end
% number of pairs
np=L*(L-1)/2;
% store log correlations between pairs with distance
V=zeros(np,2);
n=1;
for i=1:L
    for j=(i+1):L
        V(n,1)=D(i,j);
        V(n,2)=-log(C(i,j));
        n=n+1;
    end
end
% theil-sen line fit


Data.stats.C=C;

% PART 3 thouless conductances

% for each site try with H1=sigma_x_i
H1

% calculate E2_n = E_n + H1_nn
E2_d=E_d+diag(U'*H1*U);


end

