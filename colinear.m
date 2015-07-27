function [ind, p] = colinear(x, y, varargin)
%COLINEAR Find colinear points in a set of x/y coordinates
%
% [ind, p] = colinear(x, y, p1, v1, ...)
%
% Input variables:
%
%   x:  x coordinates, any size array
%
%   y:  y coordinates, same size as x
%
% Optional input arguments:
%
%   npt:    minimum number of points in a colinear set [4]
%
%   round:  rounding tolerance for angle between points to be considered
%           the same (degrees) [0.001]
%
%   thlim:  function to limit the angles considered for a set; takes as
%           input an array of theta values (in degrees) and return a
%           logical array, true to keep and false to ignore
%
% Opyright 2015 Kelly Kearney   


Opt.npt = 4;
Opt.round = 1e-5;
Opt.thlim = @(x) true(size(x));

Opt = parsepv(Opt,varargin);

nx = numel(x);

x = x(:);
y = y(:);


% lnidx = zeros(nx);
% count = 0;

sets = cell(0);
lnadj =zeros(nx);

count = 0;
for ii = 1:nx
    dx = x - x(ii);
    dy = y - y(ii);
    
    th = atand(dy./dx);
    th(th == -90) = 90;
%     th(isnan(th)) = 100;
    isn = isnan(th);
    th = th(~isn);
    thidx = find(~isn);
    
    [thunq, ith] = findclump(th, Opt.round);
    
    nper = histc(ith, 1:max(ith));
    isgrp = nper >= Opt.npt & Opt.thlim(thunq);
    
    thunqidx = find(nper > Opt.npt);

    for it = 1:length(thunqidx)
        settmp = [thidx(ith == thunqidx(it)); ii];
        [r,c] = meshgrid(settmp);
        iseq = r == c;
        idx = sub2ind([nx nx], r(~iseq), c(~iseq));
        if lnadj(idx(1)) == 0
            count = count + 1;
            lnadj(idx) = count;
            
        end 
    end
end

ind = cell(count,1);
for ii = 1:count
    [r,c] = find(lnadj == ii);
    ind{ii} = unique([r;c]);
    p(ii,:) = polyfit(x(ind{ii}), y(ind{ii}), 1);
end

        
%     
%     [srt, isrt] = sort(th);
%     
%     dth = diff(srt) == 0;
%     tmp1 = [true; dth];
%     tmp2 = [dth; true];
%     
%     se = [find(tmp1) find(tmp2)];
%     npt = se(:,2) - se(:,1) + 1;
%     se = se(npt >= Opt.npt,:);
%     
%     for iln = 1:size(se,1)
%         tmp = sort([isrt(se(iln,1):se(iln,2)); ii]);
%         
%         rc = nchoosek([tmp; ii],2);
%         rc = [rc; fliplr(rc)];
%         tmpidx = sub2ind([nx nx], rc(:,1), rc(:,2));
%         
% %         tmpidx = sort(sub2ind([nx nx], [ones(size(tmp))*ii; tmp], [tmp; ones(size(tmp))*ii]));
%         if ~lnidx(tmpidx(1))
%             count = count + 1;
%             lnidx(tmpidx) = count;
%         end
%     end
%     
% end
% 
% ind = cell(count,1);
% my = nan(count,2);
% W = warning('off');
% for ii = 1:count
%     
%     [r,c] = ind2sub([nx nx], find(lnidx == ii));
%     ind{ii} = unique([r;c]);
%     
%     
%     p(ii,:) = polyfit(x(ind{ii}), y(ind{ii}), 1);
%     
% end
% warning(W);
% 
% 


% for ii = 1:nx
%     dx = x - x(ii);
%     dy = y - y(ii);
% 
%     m(ii,:) = dy./dx;
%     
%     yint(ii,:) = y(ii) - m(ii,:).*x(ii);
%     xint(ii,:) = x(ii) + y(ii)./m(ii,:);
% end
% th = atand(m);
% th(th == -90) = 90;
% 
% isvert = th == 90; % special case, vertical
% int = yint;
% int(isvert) = xint(isvert);
% 
% [unqti, ia, ic] = unique([th(:) int(:)], 'rows');
% 
% % [unqmy0, ia, ic] = unique([m(:) y0(:)], 'rows');
% 
% nper = histc(ic, 1:max(ic));
% 
% idx = find(nper > Opt.npt);
% 
% % my = unqmy0(idx,:);
% 
% ind = cell(length(idx),1);
% my = nan(length(idx),1);
% for ii = 1:length(idx)
%     tmp = find(idx(ii) == ic);
%     [r,c] = ind2sub([nx nx], tmp);
%     ind{ii} = unique([r;c]);
%     if all(isvert(tmp))
%         my(ii,1) = Inf;
%         my(ii,2) = mean(xint(tmp));
%     else
%         my(ii,1) = mean(m(tmp));
%         my(ii,2) = mean(yint(tmp));
%     end
% end





