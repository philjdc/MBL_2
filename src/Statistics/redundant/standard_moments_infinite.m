function [mu,c,s] = standard_moments_infinite(X,M)
% calculates the moments of the data Z=(z_1,...,z_m) which is the sample X
% shift-rescaled to have mean = 0 and std.dev = 1.
% mu(n) -   nth moment
% c     -   centre (mean) of the raw distribution
% s     -   width (unbiased est of std.dev) of the raw distribution

% list data
X=X(:);
% get the center
c=mean(X);
% get the width
s=std(X,1);

% shift and rescale data to interval [-1,1]
Z=(X-c)/s;

mu=zeros(1,M);
Zm=ones(size(Z));
for m=1:M
    % get data rasied to nth power
    Zm=Zm.*Z;
    % calculate nth moment
    mu(m)=mean(Zm);
    % these estimates are biased? need to think about this a bit?
end

end

