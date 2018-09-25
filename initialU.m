function [x] =initialU(node,elem,u,ux,uy)


[~,edge,~] = dofP2(elem);
NV = size(node,1);
NT = size(elem,1);
NE = size(edge,1);
Ndof = NV + NE*5+NT;
x = zeros(Ndof,1);
x(1:NV) = u(node);
midnode = (node(edge(:,1), :) + node(edge(:,2),:))/2;
x(NV+1:NE+NV) = u(midnode );
x(NV+NE+1:4:NV+NE*5) = ux(node(edge(:,1), :));
x(NV+NE+2:4:NV+NE*5) = ux(node(edge(:,2), :));
x(NV+NE+3:4:NV+NE*5) = uy(node(edge(:,1), :));
x(NV+NE+4:4:NV+NE*5) = uy(node(edge(:,2), :));