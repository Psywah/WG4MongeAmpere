function F = WG4MongeAmpere_Input(x, elem, node,h,f,ug,ux,uy)
%% usage

% h=1;                                      % mesh size
% [node,elem] = squaremesh([0,1,0,1],h);    % mesh information
% [elem,idx,area] = fixorder(node,elem);
% [~,edge,bdDof] = dofP2(elem);
% Ndof = size(node,1)+size(edge,1)*5 +size(elem,1) - numel(bdDof) - 2*sum(float(bdDof>size(node,1)));  % number of unknows
% f = @(coord) sin(pi * coord(:,1)).* sin(pi*coord(:,2));  %source function

% x = rand(Ndof,1);
% F = WG4MongeAmpere(x, elem, node,f,h) ; 


%% generate dofs

[elem2P2dof, elem2edgeDof, elem2elemDof, edge, bdDof] = dofmap(node, elem);

%% x= recoverX(x,node,elem,edge,bdDof,ug,ux,uy);
Ndof = size(node,1)+size(edge,1)*5 +size(elem,1);
idx = 1:Ndof;
bdNodeIdx = bdDof(bdDof<=size(node,1));
bdEdgeIdx = bdDof(bdDof>size(node,1)) -size(node,1);
tmpIdx = [bdDof; bdEdgeIdx*4-3 +  size(node,1)+size(edge,1);...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+1;...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+2;...
                bdEdgeIdx*4-3 +  size(node,1)+size(edge,1)+3];
idx(tmpIdx) = [];
tmpx=sym(zeros(Ndof,1));
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
x = tmpx;




%% compute weak hessian for each basis
wh = zeros(size(elem,1), 4, 12 );
I = eye(2);
area = simplexvolume(node,elem);
tmp = [2,3; 3,1; 1,2];
for ei =1:3
    N = node(elem(:,tmp(ei,2)), :) - node(elem(:,tmp(ei,1)), :);
    N = [N(:,2), -N(:,1)];
    for bi=1:2
        for vi=1:2 
            wh(:,1,ei*4-4+bi*2-2+vi) =.5 *I(:,bi)'*I(:,1)* N(:,1) ./area;
            wh(:,2,ei*4-4+bi*2-2+vi) = .5 *I(:,bi)'*I(:,1)* N(:,2) ./area;
            wh(:,3,ei*4-4+bi*2-2+vi) = .5 *I(:,bi)'*I(:,2)* N(:,1) ./area;
            wh(:,4,ei*4-4+bi*2-2+vi) = .5 *I(:,bi)'*I(:,2)* N(:,2) ./area;
        end
    end
end



%%



ii = elem2elemDof;


whx = sym(zeros(size(elem,1),4));
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
[Dlambda,~,~] = gradbasis(node,elem);
tmp = [2,3;3,1;1,2];
for ei=1:3
    lambda1 = zeros(size(lambda,1), 3);
    lambda1(:,tmp(ei,:)) = lambda;
    N = node(elem(:,tmp(ei,2)),:) -node(elem(:,tmp(ei,1)),:) ;
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
%F = accumarray(ii,vv);
F=sym(zeros(size(unique(ii))));
for i=1:size(F,1)
    ind=find(ii==i);
    F(i)=sum(vv(ind));
end

% F(bdDof)=[];
% 
 bdEdgeIdx = bdDof(bdDof>size(node,1)) -size(node,1);
N = node(edge(bdEdgeIdx, 2), :) - node(edge(bdEdgeIdx, 1), :);
N = [N(:,2), -N(:,1)];
tmpIdx =bdEdgeIdx*4-3 +  size(node,1)+size(edge,1);
F(tmpIdx) = F(tmpIdx).*N(:,1) +     F(tmpIdx+2).*N(:,2);
F(tmpIdx+1) = F(tmpIdx+1).*N(:,1) + F(tmpIdx+3).*N(:,2);
F([bdDof;tmpIdx+2; tmpIdx+3] )=[];



end