function P = basis_function_coeffs2PDF(X,A,c,s,family)

% what order poly
M=length(A)-1;

% how many data points are there
N=length(X);

% centralise the points
X=(X-c)/s;

% get polynomial coefficients and weight functions
switch lower(family)
    case {'chebt','cheb1'}
        h=chebyshev1_poly_mat(M,true);
        W=ones(size(X));
    case 'chebf'
        h=chebyshev1_poly_mat(M,true);
        W=(1-X.^2).^(-1/4);
    case {'chebu','cheb2'}
        h=chebyshev2_poly_mat(M,true);
        W=ones(size(X));
    case 'legep'
        h=legendre_poly_mat(M,true);
        W=ones(size(X));
    case 'hermf'
        h=hermite_poly_mat(M,true);
        W=exp(-X.^2/2);
    case 'hermh'
        h=hermite_poly_mat(M,true);
        W=ones(size(X));
    otherwise
        error(['Error: did not recognise basis function "',family,'"']);
end

% constructe the vandermode matrix V_{nm} = w_n x_n^m
% - n in [1,N]
% - m in [0,M]
V=zeros(M+1,N);
V(1,:)=W;
for m=1:M
    V(m+1,:)=V(m,:).*X;
end

% evaluate distribution at these points
P=(1/s)*(A(:)')*h*V;

end

