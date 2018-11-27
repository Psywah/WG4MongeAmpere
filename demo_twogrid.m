
clear
%% define mesh size
 h=1; Nbisect = 5;
 h = h/2^(Nbisect/2);
%  f = @(coord) 4*ones(size(coord(:,1)));
% ug = @(coord) 4*zeros(size(coord(:,1)));
% ux = @(coord) 4*zeros(size(coord(:,1)));
% uy = @(coord) 4*zeros(size(coord(:,1)));

% u = x^2+ y^2 -1
 f = @(coord) 4*ones(size(coord(:,1)));
ug = @(coord) coord(:,1).^2 +coord(:,2).^2  -1;
ux = @(coord) 2*coord(:,1);
uy = @(coord) 2*coord(:,2);


%% generate mesh
% 4 elems
 noden = [0,0;1,0;1,1;0,1;0.5,0.5];
 elemn = [1,2,5;2,3,5;3,4,5;4,1,5];
 [elemn,~,~] = fixorder(noden,elemn);
 elemn = label(noden,elemn);
 for i = 1:Nbisect
    node = noden;elem=elemn;
    [noden,elemn,~,HB,tree] = bisect(node,elem,'all');
 end
%  subplot(1,2,1); showmesh(node,elem);
%  findnode(node); findelem(node,elem);% findedge(node,edge);
%  subplot(1,2,2); showmesh(noden,elemn);
%  findnode(noden); findelem(noden,elemn);% findedge(noden,edgen); 
 
 %% number of unknows
 [~,edge,bdDof] = dofP2(elem);
   [~,edgen,bdDofn] = dofP2(elemn);
 NdofU0 = size(noden,1)+size(edgen,1) - numel(bdDofn)
 %tmp = load(['sol_', int2str(2*2^Nbisect), '_ele_h1.mat']);
  tmp = load(['./data/ele', int2str(2*2^Nbisect), '_DBC.mat']);
 NinitSol= size(tmp.x,2)
 for i = 1:NinitSol
    xc(:,i)= recoverX(tmp.x(:,i),node,elem,edge,bdDof,ug,ux,uy);
 end
 xf  = WGinterpolate(xc,HB,tree,node,elem,noden,elemn);

 for i = 1:NinitSol
    x0(:,i) = BCtoX(xf(:,i),noden,elemn,edgen,bdDofn);
 end
 
 
%% sourse function



       
 fun = @(x) WG4MongeAmpere(x, elemn, noden,h,f,ug,ux,uy) ;
 
%options = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
%    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);
options = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
    'Display','final-detailed','MaxIter',100,'MaxFunEvals',1000000);
%options = optimoptions('fsolve','Algorithm','trust-region-reflective',...
%    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);
%options = optimoptions('fsolve','Algorithm','trust-region-dogleg',...
%    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);

allx=[];
for i = 1:NinitSol
    [x, F] = fsolve(fun,x0(:,i),options);
    fprintf('solution %d, resid %e\n',i,norm(F));
    if norm(F)<1e-6
        allx=[allx,x];
    end
end
tol = 1e-6;
 x = delRept(allx, tol);
 save(['./data/ele', int2str(4*2^Nbisect), '_DBC.mat'],'allx','x');
 

%%
%x = recoverX(x,node,elem,edge,bdDof);
%error_u0 = getL2error(node,elem,u,x(1:size(node,1)+size(edge,1)));
%e_u0=[e_u0 error_u0]

