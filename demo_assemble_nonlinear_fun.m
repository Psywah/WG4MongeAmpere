
%clear
%% define mesh size
 h=1; 

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

node = [0,0;1,0;1,1;0,1;0.5,0.5];  % 4 elems
elem = [1,2,5;2,3,5;3,4,5;4,1,5];

[elem,~,~] = fixorder(node,elem);

%% number of unknows
[~,edge,bdDof] = dofP2(elem);
NdofTotal = size(node,1)+size(edge,1)*5 +size(elem,1);
Ndof = NdofTotal - numel(bdDof)- 2*sum(double(bdDof>size(node,1)));
NdofU0 = size(node,1)+size(edge,1) - numel(bdDof);
disp('mesh size: ' + num2str(h));
disp(['#dof, #U0, #lambda:  ', int2str(Ndof), ', ', int2str(NdofU0), ', ', int2str(size(elem,1)) ]);



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


