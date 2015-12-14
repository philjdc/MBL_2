function [ma,mb,da,db]=theil_sen_line_fit(x,y)
% finds the gradient of a data set using the theil-sen estimator

% number of data point
n=length(x);

% number of pairs
%m=0.5*n*(n-1);

x=x(:);
y=y(:);

X=repmat(y,[1,n]);
Y=repmat(x,[1,n]);

% list of dx and dy values
dy=X-X';
dx=Y-Y';

dx=dx(:);
dy=dy(:);

% remove dx=0
I=(dx==0);
dx(I)=[];
dy(I)=[];

% list of gradients
b=dy./dx;
% find the median
mb=median(b);
% get mad
db=median(abs(b-mb));

% list of intercepts
a=y-mb*x;

% get median
ma=median(a);
% get mad
da=median(abs(a-ma));


% plot(x,y,':.');
% hold on
% xV=sort(unique(x));
% plot(xV,a+b*xV,'k')


