function [ind, p] = colinear(x, y, varargin)
%COLINEAR Find colinear points in a set of x/y coordinates
%
% [ind, p] = colinear(x, y, p1, v1, ...)
%
% In a set of scattered points, this function picks out those that fall
% along a line.
%
% Input variables:
%
%   x:      x coordinates, any size array
%
%   y:      y coordinates, same size as x
%
% Optional input arguments:
%
%   npt:    minimum number of points in a colinear set [4]
%
%   round:  rounding tolerance for angle between points to be considered
%           the same (degrees) [1e-5]
%
%   thlim:  function to limit the angles considered for a set; takes as
%           input an array of theta values (in degrees) and return a
%           logical array, true to keep and false to ignore.  For example,
%           to find only lines oriented +/-5 degrees from horizontal, use
%           @(x) abs(x) < 5.
%
% Output variables:
%
%   ind:    n x 1 cell array, where each cell holds the indices of points
%           that form a line
%
%   p:      n x 2 array, coefficients of the best-fit polynomial for each
%           set of colinear points 
%
% Copyright 2015 Kelly Kearney   


Opt.npt = 4;
Opt.round = 1e-5;
Opt.thlim = @(x) true(size(x));

Opt = parsepv(Opt,varargin);

nx = numel(x);

x = x(:);
y = y(:);

ind = cell(0);
% lnadj = zeros(nx);

% count = 0;
for ii = 1:nx
    dx = x - x(ii);
    dy = y - y(ii);
    
    th = atand(dy./dx);
    th(th == -90) = 90;
    isn = isnan(th);
    th = th(~isn);
    thidx = find(~isn);
    
    [thunq, ith] = findclump(th, Opt.round);
    
    nper = histc(ith, 1:max(ith));
    isgrp = nper >= Opt.npt & Opt.thlim(thunq);
    
    thunqidx = find(isgrp);

    for it = 1:length(thunqidx)
        settmp = sort([thidx(ith == thunqidx(it)); ii]);

        isin = cellfun(@(x) isequal(settmp, x), ind);
        if isempty(ind) || ~any(isin) 
            ind = [ind; {settmp}];
        end
        
%         [r,c] = meshgrid(settmp);
%         iseq = r == c;
%         
%         r = r(~iseq);
%         c = c(~iseq);
%         idx = sub2ind([nx nx], r(~iseq), c(~iseq));
%         if ~all(lnadj(idx)) %lnadj(idx(1)) == 0
%             count = count + 1;
%             lnadj(idx) = count;
%         end 
    end
end

nset = length(ind);
p = zeros(nset,2);

for ii = 1:nset
    p(ii,:) = polyfit(x(ind{ii}), y(ind{ii}), 1);
end



 





