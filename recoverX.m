function [fullx] = recoverX(x,node,elem,edge,bdDof)


Ndof = size(node,1)+size(edge,1)*5 +size(elem,1);
idx = 1:Ndof;
bdEdgeIdx = bdDof(bdDof>size(node,1)) -size(node,1);
tmpIdx = [bdDof; bdEdgeIdx*4-3 +  size(node,1)+size(edge,1);...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+1;...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+2;...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+3];
idx(tmpIdx) = [];
tmpx=zeros(Ndof,1);
tmpx(idx) = [x(1:end-numel(bdEdgeIdx)*2-size(elem,1));x(end-size(elem,1)+1:end)];
N = node(edge(bdEdgeIdx, 2), :) - node(edge(bdEdgeIdx, 1), :);
N = [N(:,2), -N(:,1)];
tmpx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)) = x(end-numel(bdEdgeIdx)*2-size(elem,1)+1:2:end-size(elem,1)) .*N(:,1);
tmpx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+1) = x(end-numel(bdEdgeIdx)*2-size(elem,1)+2:2:end-size(elem,1)) .*N(:,1);
tmpx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+2) = x(end-numel(bdEdgeIdx)*2-size(elem,1)+1:2:end-size(elem,1)) .*N(:,2);
tmpx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+3) = x(end-numel(bdEdgeIdx)*2-size(elem,1)+2:2:end-size(elem,1)) .*N(:,2);
fullx = tmpx;


