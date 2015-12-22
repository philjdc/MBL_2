function h=chebyshev2_poly_mat(n,normalise)
% gives a matrix of coefficients for chebyshevU polynomial
% H_n(x) = sum_{m=0}^n h(n+1,m+1) x^m
% if normalise is true these are divded by sqrt(pi/2) for n>0 and sqrt(pi)
% for n=0

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
        h(m+2,:)=[0,2*h(m+1,1:end-1)]-h(m,:);
        % increment next row
        m=m+1;
    end
end


% multiply each polynomial by 1/sqrt(pi/2) for n>0 and 1/sqrt(pi) for n=0
if normalise
    h(1,:)=h(1,:)/sqrt(pi);
    for m=1:n
        h(m+1,:)=h(m+1,:)/sqrt(pi/2);
    end
end