
clear
%% define mesh size
 h=1; Nbisect = 1;
 h = h/2^(Nbisect/2);
%  f = @(coord) 4*ones(size(coord(:,1)));
% ug = @(coord) 4*zeros(size(coord(:,1)));
% ux = @(coord) 4*zeros(size(coord(:,1)));
% uy = @(coord) 4*zeros(size(coord(:,1)));
% 
% % u = x^2+ y^2 -1
%  f = @(coord) 4*ones(size(coord(:,1)));
% u = @(coord) coord(:,1).^2 +coord(:,2).^2  -1;
% ux = @(coord) 2*coord(:,1);
% uy = @(coord) 2*coord(:,2);

% % u = 1/2*x^2+x y^2/4 
%  f = @(coord) .5*ones(size(coord(:,1)));
% u = @(coord) coord(:,1).^2/2 +coord(:,1) + coord(:,2).^2/4;
% ux = @(coord) coord(:,1)+ones(size(coord(:,1)));
% uy = @(coord) coord(:,2)/2;



% u = x^4+y^4
f = @(coord) 144*coord(:,1).^2.*coord(:,2).^2;
u = @(coord) coord(:,1).^4 +coord(:,2).^4 ;
ux = @(coord) 4*coord(:,1).^3;
uy = @(coord) 4*coord(:,2).^3;
% 
% % maping square to square
% q = @(x) (-x.^2/(8*pi) + 1/(256*pi^3) +1/(32*pi)).*cos(8*pi*x) + x.*sin(8*pi*x)/(32*pi^2);
% dq = @(x) (x.^2-.25).*sin(8*pi*x);
% ddq = @(x) 2*x.*sin(8*pi*x) + 8*pi*(x.^2-.25).*cos(8*pi*x);
% f = @(coord) ones(size(coord(:,1))) +4*(ddq(coord(:,1)).*q(coord(:,2)) +ddq(coord(:,2)).*q(coord(:,1)) )...
%                     +16*(ddq(coord(:,1)).*q(coord(:,2)).*ddq(coord(:,2)).*q(coord(:,1)) -...
%                           dq(coord(:,1)).^2.*dq(coord(:,2)).^2 );
% u = @(coord) sum(coord.^2,2)/2 + 4*q(coord(:,1)).*q(coord(:,2));
% ux = @(coord) coord(:,1) + 4*dq(coord(:,1)).*q(coord(:,2));
% uy = @(coord) coord(:,2) + 4*dq(coord(:,2)).*q(coord(:,1));


%% generate mesh
% 4 elems
 noden = [0,0;1,0;1,1;0,1;0.5,0.5];
 %noden = noden - 0.5;
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
 %tmp = load(['./data/ele', int2str(2*2^Nbisect), '_DBC.mat']);
 %tmp = load(['./data/ele', int2str(2*2^Nbisect), '_quadratic.mat']);
 tmp = load(['./data/ele', int2str(2*2^Nbisect), '_quartic.mat']);
 %tmp = load(['./data/ele', int2str(2*2^Nbisect), '_square.mat']);
 NinitSol= size(tmp.x,2)
 for i = 1:NinitSol
    xc(:,i)= recoverX(tmp.x(:,i),node,elem,edge,bdDof,u,ux,uy);
 end
 xf  = WGinterpolate(xc,HB,tree,node,elem,noden,elemn);

 for i = 1:NinitSol
    x0(:,i) = BCtoX(xf(:,i),noden,elemn,edgen,bdDofn);
 end
 
 
%% sourse function



       
 fun = @(x) WG4MongeAmpere(x, elemn, noden,h,f,u,ux,uy) ;
 
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
size(x)
%save(['./data/ele', int2str(4*2^Nbisect), '_DBC.mat'],'allx','x');
%save(['./data/ele', int2str(4*2^Nbisect), '_quadratic.mat'],'allx','x');
save(['./data/ele', int2str(4*2^Nbisect), '_quartic.mat'],'allx','x');
%save(['./data/ele', int2str(4*2^Nbisect), '_square.mat'],'allx','x');
 

%%
%x = recoverX(x,node,elem,edge,bdDof);
%error_u0 = getL2error(node,elem,u,x(1:size(node,1)+size(edge,1)));
%e_u0=[e_u0 error_u0]

