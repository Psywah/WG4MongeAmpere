function x = delRept(x, tol)

i = 1;
while i~= size(x,2)
    j = i+1;
    while j <= size(x,2)
        if norm(x(:,i)-x(:,j)) <tol
            x(:,j) =[];
        else
            j = j+1;
        end
    end
    i=i+1;
end
