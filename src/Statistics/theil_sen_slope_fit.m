function [mb,db]=theil_sen_slope_fit(x,y)
% finds the gradient of a data set using the theil-sen estimator, returns
% an estimation of the distribution width as well

% number of data point
n=length(x);

% number of pairs
m=0.5*n*(n-1);

dy=repmat(y(:),[1,n])-repmat(y(:),[1,n])';
dx=repmat(x(:),[1,n])-repmat(x(:),[1,n])';

dx=dx(:);
dy=dy(:);

% remove dx=0
I=(dx==0);
dx(I)=[];
dy(I)=[];
% convert to gradients
b=dy./dx;
% find the median
mb=median(b);
% find the median absolute deviation
db=median(abs(b-mb));


