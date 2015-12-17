function [H,n,p,m,I,Z,D] = hamiltonian_east_model(J,t,h,dh,L,left_flip,right_flip,boundary_cyclic)
% makes facilitated spin hamiltonian. Returns spin algebra operators as
% well

[n,p,m,I,Z]=manybody_algebra(L,false);

h_min=h-dh/2;
h_max=h+dh/2;

% on site field
h=unifrnd(h_min,h_max,L,1);

% initialise 
H=sparse(2^L,2^L);

% add on site fields
for i=1:L
    if abs(h(i)>0)
        H=H+h(i)*n{i};
    end
end

if (abs(J)>0)
    % add zz coupling
    for i=1:(L-1)
        H=H+J*(n{i}*n{i+1});
    end
    if boundary_cyclic
        H=H+J*(n{L}*n{1});
    end
end



if (abs(t)>0)
    % add flipping due to paticle to the left
    if left_flip
        for i=1:(L-1)
            H=H-t*n{i}*(p{i+1}+m{i+1});
        end
        if boundary_cyclic
            H=H-t*n{L}*(p{1}+m{1});
        end
    end

    % add flipping due to paticle to the right
    if right_flip
        for i=1:(L-1)
            H=H-t*n{i+1}*(p{i}+m{i});
        end
        if boundary_cyclic
            H=H-t*n{1}*(p{L}+m{L});
        end
    end
end

D=dist_matrix(adj_matrix_square_lattice(L,1,boundary_cyclic));
