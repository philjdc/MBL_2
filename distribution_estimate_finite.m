function P= distribution_estimate_finite(mu,X,x_lim)
% given the first n standard moments, the probability distribution is
% reconstructed and evaluated at the points X=(x_1,...,x_m)
% such that P_n = p(x_n)

% rescale the X values to the standard interval

% get the limits of the support
x_min=min(x_lim);
x_max=max(x_lim);

% get the center and width
c=sum(x_lim(1:2))/2;
s=abs(diff(x_lim(1:2))/2);

% get the center and width
c=sum(x_lim(1:2))/2;
s=abs(diff(x_lim(1:2))/2);
% check finite support
assert(s>0);
% shift and rescale data to interval [-1,1]
X=(X-c)/s;

% order of hughest moment
M=length(mu);

% get a legendre function matrix
h=legendre_poly_mat(M,true);

% include the 0th moment
mu=[1;mu(:)]';

% number of evaluation points
N=length(X);

% make evaluation points into column 
X=X(:);

% constructe the vandermode matrix V_{nm} = x_n^m
% - n in [1,N]
% - m in [0,M]
V=zeros(N,M+1);
V(:,1)=1;
for m=1:M
    V(:,m+1)=V(:,m).*X;
end

% weight function is trivial for legendre polynomials
W=ones(N,1);


P=W.*(V*(h')*h*(mu'));


end

