function h=hermite_poly_mat(n,normalise)
% gives a matrix of coefficients for physicists hermite polynomials
% H_n(x) = sum_{m=0}^n h(n+1,m+1) x^m
% if normalise is true these are divded by sqrt(sqrt(pi)n!2^n)

% by default just give legendre polys
if nargin==1
    normalise=false;
end

% initialise
h=zeros(n+1);

% solve by recursion
if n==0
    h=1;
else
    % initialise coefficients
    h(1,1)=1;
    h(2,2)=2;
    % initialise index last finished row
    m=1;
    while (m<n)
        h(m+2,:)=[0,2*h(m+1,1:end-1)]-2*m*h(m,:);
        % increment next row
        m=m+1;
    end
end

% multiply each polynomial by
if normalise
    for m=0:n
        h(m+1,:)=exp(-log(pi)/4-m*log(2)/2-gammaln(m+1)/2)*h(m+1,:);
    end
end