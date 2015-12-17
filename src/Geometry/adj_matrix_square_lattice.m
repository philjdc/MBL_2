function A = adj_matrix_square_lattice(L,D,periodic)
% makes an adjacency matrix for a periodic square lattice of length L in D dimensions

% number of sites
N=L^D;

A=zeros(N);

% for each site
for n=1:N
    % convert it to a coordinate vecotr
    m=n-1;
    x=zeros(1,D);
    for d=1:D
        x(d)=mod(m,L);
        m=(m-x(d))/L;
    end
    x=x+1;
    
    % for each dim
    for d=1:D
        % get two neighbours
        xL=x;
        xR=x;
        % it has two neighbours, one left and one right
        xL(d)=mod((x(d)-1)-1,L)+1;
        %disp(x-xL);
        xR(d)=mod((x(d)+1)-1,L)+1;
        %disp(x-xR);
        % convert these coordinates back to a number
        nL=1+dot((xL-1),repmat(L,1,D).^(0:D-1));
        nR=1+dot((xR-1),repmat(L,1,D).^(0:D-1));
        % put in the adjacency matrix
        if periodic
            A(n,nL)=1;
            A(nL,n)=1;
            A(n,nR)=1;
            A(nR,n)=1;
        else
            if xL(d)==(x(d)-1)
                A(n,nL)=1;
                A(nL,n)=1;
            end
            if xR(d)==(x(d)+1)
                A(n,nR)=1;
                A(nR,n)=1; 
            end
        end
    end
end



end

