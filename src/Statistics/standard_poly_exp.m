function [E_poly,c,s] = standard_poly_exp(X,W,poly_mat,x_lim)
% returns estimates of polynomials from a data sample
% standardises data
%   - if x_lim = [x_min,x_max] is specified this maps that interval to
%     [-1,1]
%   - if x_lim is not specified the data is shift-scaled to have mean = 0
%     and std.dev = 0.
% c,s give the shift and scaling used respectively
% E_poly(n) = mean( poly_n(x_i) w_i ) 
%       where poly_n(x) = sum_{m=0} poly(m+1) x^m

% standardise the data
if nargin==3
    % if infinite support shift scale by the mean std.dev
    c=mean(X);
    s=std(X,1);
else
    % else there is specified finite support - map to [-1,1]
    c=mean(x_lim)/2;
    s=abs(diff(x_lim))/2;
end
% shift scale
X=(X-c)/s;
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
%disp(size(poly_mat));
%disp(size(V));
E_poly=mean((poly_mat*V),2);

end

