function [x] = BCtoX(fullx,node,elem,edge,bdDof)




Ndof = size(node,1)+size(edge,1)*5 +size(elem,1);
idx = 1:Ndof;
bdEdgeIdx = bdDof(bdDof>size(node,1)) -size(node,1);
tmpIdx = [bdDof; bdEdgeIdx*4-3 +  size(node,1)+size(edge,1);...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+1;...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+2;...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+3];
idx(tmpIdx) = [];
Ndofx = Ndof -numel(bdDof) - 2*size(bdEdgeIdx,1);
x=zeros(Ndofx,1);
x(1:end -numel(bdEdgeIdx)*2-size(elem,1) ) =fullx(idx(1:end - size(elem,1)));
x(end - size(elem,1) +1 : end ) = fullx(end - size(elem,1) +1 : end);

T = node(edge(bdEdgeIdx, 2), :) - node(edge(bdEdgeIdx, 1), :);
T = diag(1./sqrt(sum(T.^2,2)))*T;
N = [T(:,2), -T(:,1)];
length = sum(N.^2,2);
x(end-numel(bdEdgeIdx)*2-size(elem,1)+1:2:end-size(elem,1)) =...
        (fullx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)) .*N(:,1)+...
        fullx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+2).*N(:,2))./length;
x(end-numel(bdEdgeIdx)*2-size(elem,1)+2:2:end-size(elem,1)) =...
        (fullx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+1) .*N(:,1)+...
        fullx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+3).*N(:,2))./length;


