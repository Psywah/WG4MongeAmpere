function [elem2P2dof, elem2edgeDof, elem2elemDof, edge, bdDof] = dofmap(node, elem)

[elem2P2dof,edge,bdDof] = dofP2(elem);


elem2edge = elem2P2dof(:,4:6) - size(node,1);
edge2dof = uint32( reshape(1:size(edge,1)*4, 4, size(edge,1))')+ max(elem2P2dof(:));
tmp = [2,3; 3,1; 1,2];
for i = 1:3
    elem2edgeDof(:,4*i-3:4*i) = edge2dof(elem2edge(:,i),:);
    swapidx = elem(:,tmp(i,1)) > elem(:,tmp(i,2));
    elem2edgeDof(swapidx, 4*i-3:4*i) = elem2edgeDof(swapidx, [4*i-2,4*i-3,4*i,4*i-1 ]);
end
elem2elemDof = uint32(1:size(elem,1))' + max(elem2edgeDof(:));