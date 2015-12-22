function [mu,c,s]=standard_moments_finite(X,M,x_lim)
% calculates the moments of the data Z=(z_1,...,z_m) which is the sample X
% shift-rescaled to the interval [-1,1]
% mu(n) -   nth moment
% c     -   centre of the raw distribution
% s     -   width of the raw distribution

% list data
X=X(:);
% get the limits of the support
x_min=min(x_lim);
x_max=max(x_lim);

% get the center and width
c=sum(x_lim(1:2))/2;
s=abs(diff(x_lim(1:2))/2);
% check finite support
assert(s>0);
% shift and rescale data to interval [-1,1]
Z=(X-c)/s;

mu=zeros(1,M);

Zm=ones(size(Z));
for m=1:M
    % get data rasied to nth power
    Zm=Zm.*Z;
    % calculate nth moment
    mu(m)=mean(Zm);
end

end