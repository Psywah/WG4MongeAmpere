function u = WGinterpolate(oldu,HB,tree,node,elem,noden,elemn)


T = auxstructure(elemn);
edgen = T.edge;
T = auxstructure(elem);
elem2edge = T.elem2edge; edge = T.edge;


N = size(node,1); Nn = size(noden,1);
NE = size(edge, 1); NEn = size(edgen,1);
NT = size(elem,1); NTn = size(elemn, 1);

Nf = NTn + Nn + NEn*5;
u = zeros(Nf,size(oldu,2));

%% lambda element-wise interpolate
u(Nn+NEn*5 + tree(:,2)) = oldu(N+NE*5 +tree(:,1));
u(Nn+NEn*5 + tree(:,3)) = oldu(N+NE*5 +tree(:,1));

%% u_0 and u_b
u(1:N) = oldu(1:N);
% new node
[C, ia, ib] = intersect(HB(:,[2,3]), edge,'rows','stable');
u(HB(ia,1),:) = oldu(N+ib,:);

% breaking edge 1
[C, ia, ibn] = intersect(HB(:,[2,1]), edgen,'rows','stable');
u(Nn + ibn,:) = oldu(HB(ia,2),:)*3/8  -oldu(HB(ia,3),:)/8 + oldu(N+ib,:)*3/4;
u(Nn + NEn + [ibn*4 - 3;ibn*4 - 2;ibn*4 - 1;ibn*4],:) =[ oldu(N+NE+ ib*4-3,:) ;...
                                   (oldu(N+NE+ib*4-3,:)+oldu(N+NE+ib*4-2,:))/2;...
                                    oldu(N+NE+ ib*4-1,:) ;...
                                    (oldu(N+NE+ib*4-1,:)+oldu(N+NE+ib*4,:))/2];
                                    
                                        
% breaking edge 2
[C, ia, ibn] = intersect(HB(:,[3,1]), edgen,'rows','stable');
u(Nn + ibn,:) = oldu(HB(ia,3),:)*3/8  -oldu(HB(ia,2),:)/8 + oldu(N+ib,:)*3/4;
u(Nn + NEn + [ibn*4 - 3;ibn*4 - 2;ibn*4 - 1;ibn*4],:) =[ oldu(N+NE+ ib*4-2,:) ;...
                                   (oldu(N+NE+ib*4-3,:)+oldu(N+NE+ib*4-2,:))/2;...
                                    oldu(N+NE+ ib*4,:) ;...
                                    (oldu(N+NE+ib*4-1,:)+oldu(N+NE+ib*4,:))/2];

% original edge                                
[C, ia, ib] = intersect(edge, edgen,'rows','stable');
u(Nn + ib,: ) = oldu(N+ia,:);
u(Nn+NEn + [ib*4 - 3;ib*4 - 2;ib*4 - 1;ib*4],:  ) = oldu(N+NE+[ia*4-3;ia*4-2;ia*4-1;ia*4],:);

newedge = [elem(tree(:,1),1), elemn(tree(:,2),1)];
[C, ia, ibn] = intersect(newedge, edgen,'rows','stable');


u(Nn + ibn,:) = - oldu(elem(:,2),:)/8 - oldu(elem(:,3),:)/8 + ...
              oldu(N+elem2edge(:,1),:)/4 + oldu(N+elem2edge(:,2),:)/2+ oldu(N+elem2edge(:,3),:)/2;
          
lambda = [1,0,0;0,.5,.5];
[Dlambda,~,~] = gradbasis(node,elem);
for p = 1:2
    [~, Dphip] = baseP2Dim2(lambda,Dlambda,p);
    for i = 1:2
        u(Nn+NEn + ibn*4 +(i-1)*2+p,:) =  diag(Dphip(:,i,1)) * oldu(elem(:,1),:) + ...
            diag(Dphip(:,i,2))* oldu(elem(:,2),:) + ...
            diag(Dphip(:,i,3)) * oldu(elem(:,3),:) + ...
            diag(Dphip(:,i,4)) *  oldu(N+elem2edge(:,1),:)+ ...
            diag(Dphip(:,i,5)) *  oldu(N+elem2edge(:,2),:)+...
            diag(Dphip(:,i,6)) *  oldu(N+elem2edge(:,3),:);
    end
end




