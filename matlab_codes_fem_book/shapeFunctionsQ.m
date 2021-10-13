function [shape,naturalDerivatives] = ...
    shapeFunctionsQ(xi,eta,elemType)
% shape function and derivatives for Q4, Q8 and Q9 elements
% shape: Shape functions
% naturalDerivatives: derivatives w.r.t. xi and eta
% xi, eta: natural coordinates (-1 ... +1)

switch elemType
    case 'Q4' % Q4 element
        shape = 1/4*[(1-xi)*(1-eta); (1+xi)*(1-eta);
            (1+xi)*(1+eta); (1-xi)*(1+eta)];
        
        naturalDerivatives = 1/4*[
            -(1-eta), -(1-xi); 1-eta,   -(1+xi);
            1+eta ,    1+xi; -(1+eta),  1-xi];
        
    case 'Q8' % Q8 element
        shape = 1/4*[-(1-xi)*(1-eta)*(1+xi+eta);
            -(1+xi)*(1-eta)*(1-xi+eta);
            -(1+xi)*(1+eta)*(1-xi-eta);
            -(1-xi)*(1+eta)*(1+xi-eta);
            2*(1-xi*xi)*(1-eta);
            2*(1+xi)*(1-eta*eta);
            2*(1-xi*xi)*(1+eta);
            2*(1-xi)*(1-eta*eta)];
        
        naturalDerivatives = 1/4*[
            -(eta+2*xi)*(eta-1), -(2*eta+xi)*(xi-1);
            (eta-2*xi)*(eta-1),  (2*eta-xi)*(xi+1);
            (eta+2*xi)*(eta+1),  (2*eta+xi)*(xi+1);
            -(eta-2*xi)*(eta+1), -(xi-1)*(2*eta-xi);
            4*xi*(eta-1),    2*(xi^2-1);
            2*(1-eta^2), -4*eta*(xi+1);
            -4*xi*(eta+1), 2*(1-xi^2);
            2*(eta^2-1), 4*eta*(xi-1)];
        
    case 'Q9' % Q9 element
        shape = 1/4*[xi*eta*(xi-1)*(eta-1);
            xi*eta*(xi+1)*(eta-1);
            xi*eta*(xi+1)*(eta+1);
            xi*eta*(xi-1)*(eta+1);
            -2*eta*(xi*xi-1)*(eta-1);
            -2*xi*(xi+1)*(eta*eta-1);
            -2*eta*(xi*xi-1)*(eta+1);
            -2*xi*(xi-1)*(eta*eta-1);
            4*(xi*xi-1)*(eta*eta-1)];
        
        naturalDerivatives = 1/4*[
            eta*(2*xi-1)*(eta-1),xi*(xi-1)*(2*eta-1);
            eta*(2*xi+1)*(eta-1),xi*(xi+1)*(2*eta-1);
            eta*(2*xi+1)*(eta+1),xi*(xi+1)*(2*eta+1);
            eta*(2*xi-1)*(eta+1),xi*(xi-1)*(2*eta+1);
            -4*xi*eta*(eta-1),   -2*(xi+1)*(xi-1)*(2*eta-1);
            -2*(2*xi+1)*(eta+1)*(eta-1),-4*xi*eta*(xi+1);
            -4*xi*eta*(eta+1),   -2*(xi+1)*(xi-1)*(2*eta+1);
            -2*(2*xi-1)*(eta+1)*(eta-1),-4*xi*eta*(xi-1);
            8*xi*(eta^2-1),      8*eta*(xi^2-1)];
end

end % end function