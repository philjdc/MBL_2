function [n,p,m,I,Z]=manybody_algebra(L,check)
% returns a many body algebra of L spin half particles
% n - number operator [1,0;0,0]
% p - sigma_+ [0,1;0,0]
% m - sigma_- [0,0;1,0]

% 1 site algebra
n1=sparse([1,0;0,0]);
p1=sparse([0,1;0,0]);
m1=sparse([0,0;1,0]);
I1=sparse([1,0;0,1]);

% number of sites
n=cell(1,L);
p=cell(1,L);
m=cell(1,L);
I=I1;

% initialise first site
n{1}=n1;
p{1}=p1;
m{1}=m1;

% % initialise all sites
% for i=2:L
%     for j=1:(i-1)
%        % extend other operators to new site 
%        n{j}=kron(n{j},I1);
%        p{j}=kron(p{j},I1);
%        m{j}=kron(m{j},I1);
%     end
%     % add new operators
%     n{i}=kron(I,n1);
%     p{i}=kron(I,p1);
%     m{i}=kron(I,m1);
%     
%     % extend the identity
%     I=kron(I,I1);
% end

for i=1:L
    % sites before
    nb=i-1;
    % sites after
    na=L-i;
    % dim before
    db=2^nb;
    % dim after
    da=2^na;
    % make operators
    n{i}=kron(speye(db),kron(n1,speye(da)));
    p{i}=kron(speye(db),kron(p1,speye(da)));
    m{i}=kron(speye(db),kron(m1,speye(da)));
end

% initialise (multiplicative) identity
I=speye(2^L);

% initialise zero opertors (additive identity)
Z=sparse(2^L,2^L);

if check
    % check algebra
    for i=1:L
        % operators off site commute
       for j=1:(i-1)
           assert(norm(n{i}*n{j}-n{j}*n{i},'fro')<eps);
           assert(norm(m{i}*m{j}-m{j}*m{i},'fro')<eps);
           assert(norm(p{i}*p{j}-p{j}*p{i},'fro')<eps);
           assert(norm(n{i}*p{j}-p{j}*n{i},'fro')<eps);
           assert(norm(n{i}*m{j}-m{j}*n{i},'fro')<eps);
           assert(norm(m{i}*p{j}-p{j}*m{i},'fro')<eps);
       end
       assert(norm(p{i}*m{i}-m{i}*p{i}-2*n{i}+I)<eps);
       assert(norm(n{i}*p{i}-p{i}*n{i}-p{i})<eps);
       assert(norm(n{i}*m{i}-m{i}*n{i}+m{i})<eps);
    end

    % tests passed
    disp(['algebra checks passed']);
end

end
