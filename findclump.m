function [xc, ix] = findclump(x, tol, varargin)
%FINDCLUMP Locate values close to each other
%
% [xc, ix] = findclump(x, tol, metric)
%
% This function finds groups of data points that are all within a tolerance
% value of each other and not within that tolerance to any other points in
% the dataset.
%
% Input variables:
%
%   x:          n x m array of data points, where each row represents one
%               m-dimensional point
%
%   tol:        tolerance distance
%
% Optional input variables (passed as parameter/value pairs):
%
%   metric:     distance metric to use for interpoint distance
%               1 = city block
%               2 = Euclidean (default)
%               inf = infinity-norm
%               0 = minimum difference.  
%               See ipdm help for more details.  
%
%   keepsingle: treatment of points not classified as part of a clump:
%               true:   return as part of xc, treating each point as a
%                       one-point clump 
%               false:  don't return in xc, will have an index of NaN in ix
%                       output
%
% Output values:
%
%   xc:         mean value of all isolated clumps.  This includes all
%               single data points that were not determined to be part of a
%               clump (each is considered to in its own clump).  
%
%   ix:         index of xc corresponding to each x value.

% Copyright 2015-2017 Kelly Kearney

if isvector(x)
    x = x(:);
end

p = inputParser;
p.addParameter('metric', 2, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('keepsingle', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.parse(varargin{:});

Opt = p.Results;

D = ipdm(x, 'subset', 'maximum', 'limit',  tol, 'result', 'structure', 'metric', Opt.metric);

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
if Opt.keepsingle
    ix(isn) = (1:sum(isn)) + (count - 1);
    isn = false(size(ix));
end

xc = zeros(max(ix), size(x,2));
for ii = 1:size(x,2)
    xc(:,ii) = accumarray(ix(~isn), x(~isn,ii), [max(ix) 1], @mean);
end
