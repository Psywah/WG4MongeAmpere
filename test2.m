
clear
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
Ndof = size(node,1)+size(edge,1)*5 +size(elem,1) - numel(bdDof)- 2*sum(double(bdDof>size(node,1)));
showmesh(node,elem);
findnode(node);
findedge(node,edge);

%% sourse function
f = @(coord) sin(pi * coord(:,1)).* sin(pi*coord(:,2)); 
x = rand(Ndof,1);

% nonlinear function 
F = WG4MongeAmpere(x, elem, node,f,h) ;

% %% solve
fun = @(x) WG4MongeAmpere(x, elem, node,f,h) ;
x0 = rand(Ndof,1);
options = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
    'Display','iter','MaxIter',100,'MaxFunEvals',10000);

[x, F] = fsolve(fun,x0,options)

