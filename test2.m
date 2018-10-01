
%clear
%% define mesh size
h=1;

%% generate mesh

%[node,elem] = squaremesh([0,1,0,1],h); % 2 elems

% 4 elems
 node = [0,0;1,0;1,1;0,1;0.5,0.5];
 elem = [1,2,5;2,3,5;3,4,5;4,1,5];

[elem,~,~] = fixorder(node,elem);

%% number of unknows
[~,edge,bdDof] = dofP2(elem);
NdofTotal = size(node,1)+size(edge,1)*5 +size(elem,1);
Ndof = NdofTotal - numel(bdDof)- 2*sum(double(bdDof>size(node,1)));
 showmesh(node,elem);
 findnode(node);
 findedge(node,edge);

%% sourse function
%f = @(coord) sin(pi * coord(:,1)).* sin(pi*coord(:,2));
f = @(coord) 4*(coord(:,1) -coord(:,1).^2).*(coord(:,2) -coord(:,2).^2)...
    -(ones(size(coord,1),1) - 2*coord(:,1) ).^2 .* ...
    (ones(size(coord,1),1) - 2*coord(:,2) ).^2;
f = @(p) p(:,1).^2 .* p(:,2).^2 .*(45*p(:,1).^2.*p(:,2).^2 -60*p(:,1).^2.*p(:,2)...
    +24*p(:,1).^2 - 60*p(:,1).*p(:,2).^2 +80*p(:,1).*p(:,2) -32*p(:,1) ...
    +24*p(:,2).^2 - 32*p(:,2) +12);
% u  = (x-x^2)*(y-y^2)
u = @(coord) (coord(:,1) - coord(:,1).^2).*(coord(:,2) - coord(:,2).^2);
ux = @(coord) (ones(size(coord,1),1) - 2* coord(:,1)).*(coord(:,2) - coord(:,2).^2);
uy = @(coord) (ones(size(coord,1),1) - 2* coord(:,2)).*(coord(:,1) - coord(:,1).^2);




%% nonlinear function 
x = rand(Ndof,1);
F = WG4MongeAmpere(x, elem, node,f,h) ;

%%  solve
fun = @(x) WG4MongeAmpere(x, elem, node,f,h) ;

fullx =initialU(node,elem,u,ux,uy);
x0 = BCtoX(fullx,node,elem,edge,bdDof);

options = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);

%options = optimoptions('fsolve','Algorithm','trust-region-reflective',...
%    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);
%options = optimoptions('fsolve','Algorithm','trust-region-dogleg',...
%    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);

[x, F] = fsolve(fun,x0,options);
norm(F)
%%
x = recoverX(x,node,elem,edge,bdDof);
error_u0 = getL2error(node,elem,u,x(1:size(node,1)+size(edge,1)));
%e_u0=[e_u0 error_u0]

