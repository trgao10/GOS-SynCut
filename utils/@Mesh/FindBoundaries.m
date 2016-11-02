function [BV,BE] = FindBoundaries(obj)
% find the boundary vertices of G
% out: BV: indices of the boundary vertices
% out: BE: indices of the boundary edges
% Created by Nave Zelzer on may 22 2014.
EFA = obj.ComputeEFA();
% returns edge index for every edge in EFA, thus if we have duplicates the
% indices will have duplicates too.
[~, ~, I] = unique(EFA(:,1:2), 'rows');
% count how many of each edge
binc = histc(I,unique(I));
% if only one edge in the count we have a boundary.
BE = obj.E(:,binc==1);
% get only the boundary vertices
BV = unique(BE(:))';
end