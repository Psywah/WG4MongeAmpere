
clear
%% define mesh size
 h=1; Nbisect = 4;
 h = h/2^(Nbisect/2);

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
%  findnode(node); findelem(node,elem); findedge(node,edge);
%  subplot(1,2,2); showmesh(noden,elemn);
%  findnode(noden); findelem(noden,elemn); findedge(noden,edgen); 
 
 %% number of unknows
 [~,edge,bdDof] = dofP2(elem);
   [~,edgen,bdDofn] = dofP2(elemn);
 NdofU0 = size(noden,1)+size(edgen,1) - numel(bdDofn)
 %tmp = load(['sol_', int2str(2*2^Nbisect), '_ele_h1.mat']);
  tmp = load(['sol_', int2str(2*2^Nbisect), '_ele_NS.mat']);
 NinitSol= size(tmp.x,2);
 for i = 1:NinitSol
    xc(:,i)= recoverX(tmp.x(:,i),node,elem,edge,bdDof);
 end
 xf  = WGinterpolate(xc,HB,tree,node,elem,noden,elemn);

 for i = 1:NinitSol
    x0(:,i) = BCtoX(xf(:,i),noden,elemn,edgen,bdDofn);
 end
 
 
%% sourse function
% u  = (x-x^2)*(y-y^2)
f = @(coord) 4*(coord(:,1) -coord(:,1).^2).*(coord(:,2) -coord(:,2).^2)...
                -(ones(size(coord,1),1) - 2*coord(:,1) ).^2 .* ...
                  (ones(size(coord,1),1) - 2*coord(:,2) ).^2; 
% u  = (x^2-x^3)*(y^2-y^3)
% f = x^2*y^2*(45*x^2*y^2 - 60*x^2*y + 24*x^2 - 60*x*y^2 + 80*x*y - 32*x + 24*y^2 - 32*y + 12)
f = @(p) p(:,1).^2 .* p(:,2).^2 .*(45*p(:,1).^2.*p(:,2).^2 -60*p(:,1).^2.*p(:,2)...
           +24*p(:,1).^2 - 60*p(:,1).*p(:,2).^2 +80*p(:,1).*p(:,2) -32*p(:,1) ...
           +24*p(:,2).^2 - 32*p(:,2) +12);
       
 fun = @(x) WG4MongeAmpere(x, elemn, noden,f,h) ;
 
%options = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
%    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);
options = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
    'Display','final-detailed','MaxIter',100,'MaxFunEvals',1000000);
%options = optimoptions('fsolve','Algorithm','trust-region-reflective',...
%    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);
%options = optimoptions('fsolve','Algorithm','trust-region-dogleg',...
%    'Display','iter','MaxIter',100,'MaxFunEvals',1000000);

for i = 1:NinitSol
    [allx(:,i), F] = fsolve(fun,x0(:,i),options);
    fprintf('solution %d, resid %e\n',i,norm(F));
end
tol = 1e-6;
 x = delRept(allx, tol);
 save(['sol_' int2str(4*2^Nbisect) '_ele_NS.mat'],'allx','x');
 

%%
%x = recoverX(x,node,elem,edge,bdDof);
%error_u0 = getL2error(node,elem,u,x(1:size(node,1)+size(edge,1)));
%e_u0=[e_u0 error_u0]

