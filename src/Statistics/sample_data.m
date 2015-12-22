function Data = sample_data(X,N,x_lim)
% stores data on a probability distribution including basis function
% coefficients up to order N.

% centre
Data.C.med=median(X);
Data.C.mean=mean(X);
% dispersion
Data.D.iqr=iqr(X);
Data.D.std=std(X);
Data.D.max=max(X);
Data.D.min=min(X);
Data.D.mad=median(abs(X-median(X)));
Data.D.percentile=quantile(X,100);



if nargin==2
    Data.D.finite_support=false;
    % therefore support is infinite
    % get hermite function coeffs
    [AH,c,s]=sample2basis_function_coeffs(X,N,'hermf');
    Data.B.HF=AH;
else
    % support is finite
    Data.D.finite_support=true;
    % get chebeyshev T coeffs
    [AT,c,s]=sample2basis_function_coeffs(X,N,'chebt',x_lim);
    Data.B.CT=AT;
    % get legendre P coeffs
    [AL,c,s]=sample2basis_function_coeffs(X,N,'legep',x_lim);
    Data.B.LP=AL;
    
end
Data.B.shift=c;
Data.B.scale=s;


end

