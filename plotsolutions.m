%
clear;
disp('plot solutions');
%% define mesh size
 h=1; Nbisect = 4;
  h = h/2^(Nbisect/2);
%   f = @(coord) 4*ones(size(coord(:,1)));
% u = @(coord) 4*zeros(size(coord(:,1)));
% ux = @(coord) 4*zeros(size(coord(:,1)));
% uy = @(coord) 4*zeros(size(coord(:,1)));
 f = @(coord) 4*ones(size(coord(:,1)));
u = @(coord) coord(:,1).^2 +coord(:,2).^2  -1;
ux = @(coord) 2*coord(:,1);
uy = @(coord) 2*coord(:,2);

% maping square to square
q = @(x) (-x.^2/(8*pi) + 1/(256*pi^3) +1/(32*pi)).*cos(8*pi*x) + x.*sin(8*pi*x)/(32*pi^2);
dq = @(x) (x.^2-.25).*sin(8*pi*x);
ddq = @(x) 2*x.*sin(8*pi*x) + 8*pi*(x.^2-.25).*cos(8*pi*x);
f = @(coord) ones(size(coord(:,1))) +4*(ddq(coord(:,1)).*q(coord(:,2)) +ddq(coord(:,2)).*q(coord(:,1)) )...
                    +16*(ddq(coord(:,1)).*q(coord(:,2)).*ddq(coord(:,2)).*q(coord(:,1)) -...
                          dq(coord(:,1)).^2.*dq(coord(:,2)).^2 );
u = @(coord) sum(coord.^2,2)/2 + 4*q(coord(:,1)).*q(coord(:,2));
ux = @(coord) coord(:,1) + 4*dq(coord(:,1)).*q(coord(:,2));
uy = @(coord) coord(:,2) + 4*dq(coord(:,2)).*q(coord(:,1));

% [node,elem] = squaremesh([0,1,0,1],1/100); % 2 elems
% node = node - 0.5;
% showsolution(node,elem,u(node));



%% generate mesh
% 4 elems
 noden = [0,0;1,0;1,1;0,1;0.5,0.5];
  noden = noden - 0.5;
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
 [~,edgen,bdDofn] = dofP2(elemn);
 tmp = load(['./data/ele', int2str(4*2^Nbisect), '_square.mat']);

%  
%  for i = 1:size(tmp.allx,2)
%     allx(:,i)= recoverX(tmp.allx(:,i),node,elem,edge,bdDof,u,ux,uy);
%  end
%   N = size(node,1);
%  tx = delRept(allx(1:N,:), 1e-3);
%  
 
%  
%  for i = 1:size(tmp.allx,2)
%     allx(:,i)= recoverX(tmp.allx(:,i),noden,elemn,edgen,bdDofn,u,ux,uy);
%  end
%   N = size(noden,1);
%  tx = delRept(allx(1:N,:), 1e-3);

% for i =1:size(tmp.x,2)
%     tx(:,i)=recoverX(tmp.x(:,i),node,elem,edge,bdDof,u,ux,uy);
% end
% 
for i =1:size(tmp.x,2)
    tx(:,i)=recoverX(tmp.x(:,i),noden,elemn,edgen,bdDofn,u,ux,uy);
end

  NSol= size(tx,2);
 Nplot = ceil(sqrt(NSol));
 for i = 1 :NSol
     subplot(Nplot,Nplot, i);
     showsolution(noden,elemn,tx(1:size(noden,1),i));
     zlim([0,0.3]);
     error_u0 = getL2error(noden,elemn,u,tx(1:size(noden,1)+size(edgen,1), i))
     title(int2str(i));
 end
 fullx =initialU(noden,elemn,u,ux,uy);
 
 for i =1 : NSol
     for j = i+1:NSol
         dis(i,j) = norm(tx(1:1:size(noden,1)+size(edgen,1),i)-tx(1:1:size(noden,1)+size(edgen,1),j));
     end
 end

 for i =1:NSol
     dis2u(i)= max(abs(fullx(1:size(noden,1)+size(edgen,1))-tx(1:size(noden,1)+size(edgen,1),i)));
 end
 
 figure;
  showsolution(noden,elemn,fullx(1:size(noden,1)));
 
 %subplot(Nplot,Nplot, NSol+1);
 %title([int2str(4*2^Nbisect) 'elems',  int2str(NSol) 'sols'])
 
 %x = recoverX(x,node,elem,edge,bdDof);
%error_u0 = getL2error(noden,elemn,u,tx(1:size(noden,1)+size(edgen,1)));

 

