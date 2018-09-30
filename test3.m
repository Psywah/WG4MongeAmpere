
%clear
%% define mesh size
h=1/4;

%% generate mesh

%[node,elem] = squaremesh([0,1,0,1],h); % 2 elems

% 4 elems
 node = [0,0;1,0;1,1;0,1;0.5,0.5];
 elem = [1,2,5;2,3,5;3,4,5;4,1,5];
 
 [elem,~,~] = fixorder(node,elem);
 elem = label(node,elem);
 
 %% number of unknows
 [~,edge,bdDof] = dofP2(elem);
 subplot(1,2,1);
 showmesh(node,elem);
 findnode(node);
 findelem(node,elem);
 findedge(node,edge);
 
 
 NdofTotal = size(node,1)+size(edge,1)*5 +size(elem,1);
 Ndof = NdofTotal - numel(bdDof)- 2*sum(double(bdDof>size(node,1)));
 
 
 [noden,elemn,~,HB,tree] = bisect(node,elem,'all');
  [~,edgen,bdDofn] = dofP2(elemn);
 subplot(1,2,2);
 showmesh(noden,elemn);
 findnode(noden);
 findelem(noden,elemn);
 findedge(noden,edgen);
 
 NdofTotaln = size(noden,1)+size(edgen,1)*5 +size(elemn,1);
 Ndofn = NdofTotaln - numel(bdDofn)- 2*sum(double(bdDofn>size(noden,1)));
 


%% sourse function
%f = @(coord) sin(pi * coord(:,1)).* sin(pi*coord(:,2)); 
f = @(coord) 4*(coord(:,1) -coord(:,1).^2).*(coord(:,2) -coord(:,2).^2)...
                -(ones(size(coord,1),1) - 2*coord(:,1) ).^2 .* ...
                  (ones(size(coord,1),1) - 2*coord(:,2) ).^2; 
% u  = (x-x^2)*(y-y^2)
u = @(coord) (coord(:,1) - coord(:,1).^2).*(coord(:,2) - coord(:,2).^2);
ux = @(coord) (ones(size(coord,1),1) - 2* coord(:,1)).*(coord(:,2) - coord(:,2).^2);
uy = @(coord) (ones(size(coord,1),1) - 2* coord(:,2)).*(coord(:,1) - coord(:,1).^2);


oldx =initialU(node,elem,u,ux,uy);

newx  = WGinterpolate(oldx,HB,tree,node,elem,noden,elemn);




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

