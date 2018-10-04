function err = getL2errorEdge(node,elem,ux,uy,x)


[~,edge,~] = dofP2(elem);
NE= size(edge,1);
[lambda,weight] = quadpts1(5);

nQuad = size(lambda,1);
err = zeros(NE,1);
for p = 1:nQuad
    % evaluate uh at quadrature point
    uh = [ x(1:4:4*NE)* lambda(p,1) + x(2:4:4*NE)* lambda(p,2), ...
           x(3:4:4*NE)* lambda(p,1) + x(4:4:4*NE)* lambda(p,2)];
    coord = node(edge(:,1), :)*lambda(p,1) +  node(edge(:,2), :)*lambda(p,2);
    uexact = [ ux(coord), uy(coord)];
    
    err = err + weight(p) * sum((uh - uexact).^2,2);
end
length = sqrt(sum((node(edge(:,1),:) - node(edge(:,2),:)).^2,2) );
err = sum(length.*err);

