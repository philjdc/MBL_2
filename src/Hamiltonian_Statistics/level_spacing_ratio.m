function [mr,r] = level_spacing_ratio(E)

if ~issroted(E)
    E=sort(E);
end

% find level diffs
dE=E(2:end)-E(1:end-1);
% find ratios
r=dE(1:end-1)./dE(2:end);
% take min of rk and 1/rk so that rk in [0,1]
r=min(r,r.^(-1));
% take average level spacing ratio
mr=mean(r);

end

