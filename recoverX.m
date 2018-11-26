function [fullx] = recoverX(x,node,elem,edge,bdDof,ug,ux,uy)


Ndof = size(node,1)+size(edge,1)*5 +size(elem,1);
idx = 1:Ndof;
bdNodeIdx = bdDof(bdDof<=size(node,1));
bdEdgeIdx = bdDof(bdDof>size(node,1)) -size(node,1);
tmpIdx = [bdDof; bdEdgeIdx*4-3 +  size(node,1)+size(edge,1);...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+1;...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+2;...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+3];
idx(tmpIdx) = [];
tmpx=zeros(Ndof,1);
tmpx(idx) = [x(1:end-numel(bdEdgeIdx)*2-size(elem,1));x(end-size(elem,1)+1:end)];

bdNode = [node(bdNodeIdx,:); (node(edge(bdEdgeIdx,1),: )+node(edge(bdEdgeIdx,2),: ))/2];
tmpx(bdDof) = ug(bdNode);



T = node(edge(bdEdgeIdx, 2), :) - node(edge(bdEdgeIdx, 1), :);
T = diag(1./sqrt(sum(T.^2,2)))*T;
N = [T(:,2), -T(:,1)];
grad1 = [ux(node(edge(bdEdgeIdx,1),:)),uy(node(edge(bdEdgeIdx,1),:))];
T1 = diag(sum(grad1.*T,2))*T;
grad2 = [ux(node(edge(bdEdgeIdx,2),:)),uy(node(edge(bdEdgeIdx,2),:))];
T2 = diag(sum(grad2.*T,2))*T;


tmpx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)) = x(end-numel(bdEdgeIdx)*2-size(elem,1)+1:2:end-size(elem,1)) .*N(:,1) + T1(:,1);
tmpx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+1) = x(end-numel(bdEdgeIdx)*2-size(elem,1)+2:2:end-size(elem,1)) .*N(:,1) + T2(:,1);
tmpx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+2) = x(end-numel(bdEdgeIdx)*2-size(elem,1)+1:2:end-size(elem,1)) .*N(:,2)+T1(:,2);
tmpx(bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+3) = x(end-numel(bdEdgeIdx)*2-size(elem,1)+2:2:end-size(elem,1)) .*N(:,2) + T2(:,2);
fullx = tmpx;


