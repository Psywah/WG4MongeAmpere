%
clear;
disp('plot solutions');
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
 [~,edgen,bdDofn] = dofP2(elemn);
 tmp = load(['sol_f4_elem', int2str(4*2^Nbisect), '.mat']);
 for i = 1:size(tmp.allx,2)
    allx(:,i)= recoverX(tmp.allx(:,i),noden,elemn,edgen,bdDofn);
 end
  N = size(noden,1);
 tx = delRept(allx(1:N,:), 1e-3);

% for i =1:size(tmp.x,2)
%     tx(:,i)=recoverX(tmp.x(:,i),noden,elemn,edgen,bdDofn);
% end
  NSol= size(tx,2);
 Nplot = ceil(sqrt(NSol));
 for i = 1 :NSol
     subplot(Nplot,Nplot, i);
     showsolution(noden,elemn,tx(1:size(noden,1),i));
     title(int2str(i));
 end
 %subplot(Nplot,Nplot, NSol+1);
 %title([int2str(4*2^Nbisect) 'elems',  int2str(NSol) 'sols'])
 

