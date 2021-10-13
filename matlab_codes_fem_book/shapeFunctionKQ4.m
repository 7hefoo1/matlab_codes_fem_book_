function [shape,naturalDerivatives] = shapeFunctionKQ4(xi,eta)

% shape function and derivatives for Q4 elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi and eta
% xi, eta: natural coordinates (-1 ... +1)

shape=1/4*[ (1-xi)*(1-eta);(1+xi)*(1-eta);
    (1+xi)*(1+eta);(1-xi)*(1+eta)];

% natural derivatives order:
% [d/dx, d/dy, d^2/dx^2, d^2/dy^2, d^2/dxdy]
naturalDerivatives=...
    1/4*[-(1-eta), -(1-xi);1-eta,    -(1+xi);
    1+eta,      1+xi;-(1+eta),   1-xi];
naturalDerivatives(:,5)=[1/4;  -1/4;  1/4; -1/4];

end 

