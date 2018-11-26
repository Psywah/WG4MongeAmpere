
clear

%% generate mesh
h=1;
[node,elem] = squaremesh([0,1,0,1],h);
[elem,idx,area] = fixorder(node,elem);


%% generate dofs
[elem2P2dof,edge,bdDof] = dofP2(elem);
Ndof = size(node,1)+size(edge,1)*5 +size(elem,1);
elem2edge = elem2P2dof(:,4:6) - size(node,1);
edge2dof = uint32( reshape(1:size(edge,1)*4, 4, size(edge,1))')+ max(elem2P2dof(:));
tmp = [2,3; 3,1; 1,2];
for i = 1:3
    elem2edgeDof(:,4*i-3:4*i) = edge2dof(elem2edge(:,i),:);
    swapidx = elem(:,tmp(i,1)) > elem(:,tmp(i,2));
    elem2edgeDof(swapidx, 4*i-3:4*i) = elem2edgeDof(swapidx, [4*i-2,4*i-3,4*i,4*i-1 ]);
end
elem2elemDof = uint32(1:size(elem,1))' + max(elem2edgeDof(:));

%% compute weak hessian for each basis
wh = zeros(size(elem,1), 4, 12 );
I = eye(2);
area = simplexvolume(node,elem);
tmp = [2,3; 3,1; 1,2];
for ei =1:3
    N = node(elem(:,tmp(i,2)), :) - node(elem(:,tmp(i,1)), :);
    N = [N(:,2), -N(:,1)];
    for bi=1:2
        for vi=1:2
            wh(:,1,ei*4-4+bi*2-2+vi) = .5 *I(:,bi)'*I(:,1)* N(:,1) ./area;
            wh(:,2,ei*4-4+bi*2-2+vi) = .5 *I(:,bi)'*I(:,1)* N(:,2) ./area;
            wh(:,3,ei*4-4+bi*2-2+vi) = .5 *I(:,bi)'*I(:,2)* N(:,1) ./area;
            wh(:,4,ei*4-4+bi*2-2+vi) = .5 *I(:,bi)'*I(:,2)* N(:,2) ./area;
        end
    end
end



%%


F = zeros(Ndof,1);
f = @(coord) sin(pi * coord(:,1)).* sin(pi*coord(:,2)); 

ii = elem2elemDof;

x = rand(Ndof,1);
whx = zeros(size(elem,1),4);
for i = 1: 12
    for j =1:4
        whx(:,j) = whx(:,j) + wh(:,j,i).* x(elem2edgeDof(:,i));
    end
end
midnode = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
vv = (whx(:,1).*whx(:,4) - whx(:,2).*whx(:,3)) .* area - f(midnode).*area ;

cofWhx = whx(:,[4,3,2,1]);
cofWhx(:,[2,3]) = -cofWhx(:,[2,3]);

for ei = 1:3
    for bi =1:2
        for vi=1:2
            ii = [ii; elem2edgeDof(:,ei*4-4+bi*2-2+vi)];
            vv = [vv; sum(cofWhx.*squeeze(wh(:,:,ei*4-4+bi*2-2+vi)),2) .*x(elem2elemDof).*area];
        end
    end
end

I = eye(2);
[lambda,weight] = quadpts1(2);
[Dlambda,area,elemSign] = gradbasis(node,elem);
tmp = [2,3;3,1;1,2];
for ei=1:3
    lambda1 = zeros(size(lambda,1), 3);
    lambda1(:,tmp(ei,:)) = lambda;
    N = node(elem(:,tmp(ei,2))) -node(elem(:,tmp(ei,1))) ;
    length = sqrt(sum(N.^2,2));
    for p =1:size(weight)
        [~, Dphip] = baseP2Dim2(lambda1,Dlambda,p);
        Du = zeros(size(elem,1),2);
        for bi =1:6
            Du = Du +diag(x(elem2P2dof(:,bi)))*Dphip(:,:,bi);
        end
        for bi =1:2
            for vi = 1:2
                Du = Du - x(elem2edgeDof(:,ei*4-4+bi*2-2+vi)) *lambda(p,vi)*I(bi,:);
            end
        end     
        for bi = 1:2
            for vi = 1:2
                ii = [ ii; elem2edgeDof(:,ei*4-4+bi*2-2+vi)];
                vv =[vv; ...
                sum(-Du .*( ones(size(elem,1),1)*lambda(p,vi)*I(bi,:)),2)/h .*length *weight(p)];               
            end
        end
        for bi = 1:6
            ii = [ii; elem2P2dof(:,bi)];
            vv =[vv; ...
                sum(Du .* Dphip(:,:,bi),2)/h .*length *weight(p)];
        end       
    end
end
F = accumarray(ii,vv);




