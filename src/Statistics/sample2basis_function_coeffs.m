function [A,c,s] = sample2basis_function_coeffs(X,N,family,x_lim)
% estimates basis coefficients for various families of basis functions from
% a sample of data. The data is standardised before estimation.
% The estimators are not unbiased and so for small sample
% the negation of Bessel's correction will cause errors. If unbiased
% estimators are required - need to look up h-statstics.
% Inputs
%   X       -   sample data
%   family  -   basis function family, either 'ChebT', 'ChebU', 'LegeP' or
%               'HermH'
%   N       -   order of basis functions to go to
%   x_lim   -   optional argument, if support is finite give the bounds. 
% Outputs 
%   A       -   coefficients
%   c       -   shift
%   s       -   scaling
% Basis function families
%   ChebT   -   Chebyshev polynomial of the first kind
%   ChebF   -   Chebyshev functions of the first kind (weight function
%               split to make them into a linear basis set).
%   ChebU   -   Chebyshev polynomial of the second kind
%   LegeP   -   Legendre polynomial
%   HermH   -   Hermite polynomial
%   HermF   -   Hermite functions


% standardise the data
if nargin==3
    % if infinite support shift scale by the mean std.dev
    c=mean(X);
    s=std(X,1);
else
    % else there is specified finite support - map to [-1,1]
    c=mean(x_lim);
    s=abs(diff(x_lim))/2;
    assert(max(X)<=max(x_lim));
    assert(min(X)>=min(x_lim));
end

X=(X-c)/s;


% get polynomial coefficients and weight functions
switch lower(family)
    case {'chebt','cheb1'}
        h=chebyshev1_poly_mat(N,true);
        W=(1-X.^2).^(-1/2);
    case 'chebf'
        h=chebyshev1_poly_mat(N,true);
        W=(1-X.^2).^(-1/4);
    case {'chebu','cheb2'}
        h=chebyshev2_poly_mat(N,true);
        W=(1-X.^2).^(1/2);
    case 'legep'
        h=legendre_poly_mat(N,true);
        W=ones(size(X));
    case 'hermf'
        h=hermite_poly_mat(N,true);
        W=exp(-X.^2/2);
    case 'hermh'
        h=hermite_poly_mat(N,true);
        W=exp(-X.^2);
    otherwise
        error(['Error: did not recognise basis function "',family,'"']);
end

% get coefficients
A=poly_exp(X,W,h);


end

