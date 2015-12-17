function D = dist_matrix(A)
% given an adjacency matrix computes a distance matrix

% makes sure only ones in the adjacency matrix
A=1*(A>0);

% num sites
N=size(A,1);

% lowest impossible distance
dimp=N;

% initialise to impossibly large distance except on diag
D=(ones(N)-eye(N))*dimp;

% initialise power of A
Ad=A;
% and counter 
d=1;

while any(D(:)==dimp)
    % get power of A, set zeros to dimp and ones to d
    B=(Ad==0)*dimp+(Ad>0)*d;
    % update D
    D=min(D,B);
    % update A and d;
    Ad=Ad*A;
    d=d+1;
end
    
    
