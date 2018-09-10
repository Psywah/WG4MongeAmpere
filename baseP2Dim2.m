function [phi, Dphip] = baseP2Dim2(lambda,Dlambda,p)
    phi=[];
    Dphip = [];
    if p ==0
        phi(:,1) = lambda(:,1) .* (2*lambda(:,1) - 1);
        phi(:,2) = lambda(:,2) .* (2*lambda(:,2) - 1);
        phi(:,3) = lambda(:,3) .* (2*lambda(:,3) - 1);
        phi(:,4) = 4 * lambda(:,2) .* lambda(:,3);
        phi(:,5) = 4 * lambda(:,1) .* lambda(:,3);
        phi(:,6) = 4 * lambda(:,1) .* lambda(:,2);
    else
        Dphip = [];
        Dphip(:,:,1) = Dlambda(:,:,1)*(4*lambda(p,1)-1);
        Dphip(:,:,2) = Dlambda(:,:,2)*(4*lambda(p,2)-1);
        Dphip(:,:,3) = Dlambda(:,:,3)*(4*lambda(p,3)-1);
        Dphip(:,:,4) = 4*(Dlambda(:,:,3)*lambda(p,2) + Dlambda(:,:,2)*lambda(p,3));
        Dphip(:,:,5) = 4*(Dlambda(:,:,1)*lambda(p,3) + Dlambda(:,:,3)*lambda(p,1));
        Dphip(:,:,6) = 4*(Dlambda(:,:,1)*lambda(p,2) + Dlambda(:,:,2)*lambda(p,1));
    end
end

