function x=newton_bertini(elem, node,h,f,u,ux,uy,x0)
%x0 is the initial guess after interpolation
y=[];
NN=length(x0);
for i=1:NN
    y=[y;sym(['y', num2str(i)])];
end
fid=fopen('input','wt');

%% Configure part
fprintf(fid,'CONFIG\n');
fprintf(fid,'MPTYPE:0;\n');
fprintf(fid,'USERHOMOTOPY:1;\n');
fprintf(fid,'FINALTOL: 1e-16;\n');
fprintf(fid,'END;\n\n');

%% Inputs
fprintf(fid,'INPUT\n');
ff=['function '];
for j=1:NN
    if j==NN
        ff=[ff 'f',num2str(j)];
    else
        ff=[ff 'f',num2str(j),','];
    end
end
variable=['variable '];
for j=1:NN-1
    variable=[variable,char(y(j)),','];
end
variable=[variable,char(y(NN))];
fprintf(fid,[variable,';\n']);
fprintf(fid,[ff,';\n']);

fprintf(fid,['pathvariable s;\n']);
fprintf(fid,['parameter t;\n']);
fprintf(fid,['t=s;\n']);
%F=WG4MongeAmpere_Input(y, elem, node,f,h);
F = WG4MongeAmpere_Input(y, elem, node,h,f,u,ux,uy) ;
fun = @(x) WG4MongeAmpere(x, elem, node,h,f,u,ux,uy) ;
F0=fun(x0);
for j=1:NN
    ftmp = ['f',num2str(j), ' = ',char(vpa(F(j))),'-t*(' num2str(F0(j),'%1.16e\n') ');\n'];
    fprintf(fid,ftmp);
end
fclose(fid);
fid=fopen('start','wt');
fprintf(fid,'1\n');
fprintf(fid,'%1.16f  %1.16f\n',[x0 zeros(size(x0))]');
fclose(fid);

unix('./bertini1 input')

fid=fopen('real_solutions','r');
num=fscanf(fid,'%e',1);
x=[];
for i=1:num
    xtmp=[];
    for j=1:NN
        xtmp(j)=fscanf(fid,'%e',1);
        tmp=fscanf(fid,'%e',1);
    end
    xtmp=xtmp';
    isinx=0;
    for j=1:size(x,2)
        if norm(x(:,j)-xtmp)<1e-10
            isinx=1;
            break
        end
    end
    if ~isinx
        x=[x xtmp];
    end
end
fclose(fid);
