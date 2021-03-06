function [E_poly] = poly_exp(X,W,poly_mat)
% calculates wieghted expectation values of polynomials from a data sample
% standardises data
%   - if x_lim = [x_min,x_max] is specified this maps that interval to
%     [-1,1]
%   - if x_lim is not specified the data is shift-scaled to have mean = 0
%     and std.dev = 0.
% c,s give the shift and scaling used respectively
% calculate weighted polynomial expectation values
% E_poly(n) = mean( poly_n(x_i) w_i ) 
%       where poly_n(x) = sum_{m=0} poly(m+1) x^m

X=X(:)';
W=W(:);

% get max order of poly(n)
M=size(poly_mat,2)-1;

% amount of data
N=length(X);

% constructe the vandermode matrix V_{nm} = w_n x_n^m
% - n in [1,N]
% - m in [0,M]
V=zeros(M+1,N);
V(1,:)=W;
for m=1:M
    V(m+1,:)=V(m,:).*X;
end
[~,I]=sort(X);
M=(poly_mat*V)';
%plot(X(I),M(I,:))
E_poly=mean((poly_mat*V),2);

end

