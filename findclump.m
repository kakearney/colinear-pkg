function [xc, ix] = findclump(x, tol)
%FINDCLUMP Locate values close to each other
%
% [xc, ix] = findclump(x, tol)
%
% This function finds groups of data points that are all within a tolerance
% value of each other and not within that tolerance to any other points in
% the dataset.
%
% Input variables:
%
%   x:      array of data points
%
%   tol:    tolerance value
%
% Output values:
%
%   xc:     mean value of all isolated clumps, including single data points
%           that are not clumped together.
%
%   ix:     index of xc corresponding to each x value.

% Copyright 2015 Kelly Kearney

x = x(:);

D = ipdm(x, 'subset', 'maximum', 'limit',  tol, 'result', 'structure');

ix = nan(size(x));
count = 1;
for ii = 1:length(x)
    pot = D.columnindex(D.rowindex == ii);
    
    tf = ismember(D.columnindex, pot);
    isclump = all(ismember(D.rowindex(tf), [pot; ii]));
   
    if isclump && isnan(ix(ii))
        ix([pot; ii]) = count;
        count = count + 1;
    end
    
end

isn = isnan(ix);
ix(isn) = (1:sum(isn)) + (count - 1);

xc = accumarray(ix, x, [max(ix) 1], @mean);
