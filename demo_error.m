
%clear
%% define mesh size
N=8;
h=1/N; 

%% sourse function
%f = @(coord) sin(pi * coord(:,1)).* sin(pi*coord(:,2));

% % u  = (x-x^2)*(y-y^2)
% f = @(coord) 4*(coord(:,1) -coord(:,1).^2).*(coord(:,2) -coord(:,2).^2)...
%     -(ones(size(coord,1),1) - 2*coord(:,1) ).^2 .* ...
%     (ones(size(coord,1),1) - 2*coord(:,2) ).^2;
% u = @(coord) (coord(:,1) - coord(:,1).^2).*(coord(:,2) - coord(:,2).^2);
% ux = @(coord) (ones(size(coord,1),1) - 2* coord(:,1)).*(coord(:,2) - coord(:,2).^2);
% uy = @(coord) (ones(size(coord,1),1) - 2* coord(:,2)).*(coord(:,1) - coord(:,1).^2);

% u = x^2+ y^2 -1
 f = @(coord) 4*ones(size(coord(:,1)));
u = @(coord) coord(:,1).^2 +coord(:,2).^2  -1;
ux = @(coord) 2*coord(:,1);
uy = @(coord) 2*coord(:,2);



%% generate mesh

[node,elem] = squaremesh([0,1,0,1],h); % 2 elems

[elem,~,~] = fixorder(node,elem);

%% number of unknows
[~,edge,bdDof] = dofP2(elem);
NdofTotal = size(node,1)+size(edge,1)*5 +size(elem,1);
Ndof = NdofTotal - numel(bdDof)- 2*sum(double(bdDof>size(node,1)));
NdofU0 = size(node,1)+size(edge,1) - numel(bdDof);
disp(['mesh size: ', num2str(h)]);
disp(['#dof, #U0, #lambda:  ', int2str(Ndof), ', ', int2str(NdofU0), ', ', int2str(size(elem,1)) ]);
%  showmesh(node,elem);
%  findnode(node);
%  findedge(node,edge);


%% nonlinear function 
x = rand(Ndof,1);
F = WG4MongeAmpere(x, elem, node,h,f,u,ux,uy) ;

%%  solve
fun = @(x) WG4MongeAmpere(x, elem, node,h,f,u,ux,uy) ;

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
x = recoverX(x,node,elem,edge,bdDof,u,ux,uy);
error_u0 = getL2error(node,elem,u,x(1:size(node,1)+size(edge,1)))
error_ug = getL2errorEdge(node,elem,ux,uy,x(size(node,1)+size(edge,1)+1:end-size(elem,1 )))
% e_ug = [e_ug error_ug]
% if size(e_ug,2)>1
% order = log(e_ug(1:end-1)./e_ug(2:end)) /log(2)
% end
% e_u0=[e_u0 error_u0]
% if size(e_ug,2)>1
% order = log(e_u0(1:end-1)./e_u0(2:end)) /log(2)
% end

