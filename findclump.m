function [xc, ix] = findclump(x, tol, metric)
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
%   x:      n x m array of data points, where each row represents one
%           m-dimensional point
%
%   tol:    tolerance distance
%
%   metric: distance metric to use for interpoint distance
%           1 = city block, 2 = Euclidean (default), inf = infinity-norm, 0
%           = minimum difference.  See ipdm help for more details.  
%
% Output values:
%
%   xc:     mean value of all isolated clumps, including single data points
%           that are not clumped together.
%
%   ix:     index of xc corresponding to each x value.

% Copyright 2015 Kelly Kearney

if isvector(x)
    x = x(:);
end

if nargin < 3
    metric = 2;
end

D = ipdm(x, 'subset', 'maximum', 'limit',  tol, 'result', 'structure', 'metric', metric);

ix = nan(size(x,1),1);
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

xc = zeros(max(ix), size(x,2));
for ii = 1:size(x,2)
    xc(:,ii) = accumarray(ix, x(:,ii), [max(ix) 1], @mean);
end
