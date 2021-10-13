function [shape,naturalDerivatives] = shapeFunctionK12(xi,eta)
% shape function and derivatives for not conforming elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi and eta
% xi, eta: natural coordinates (-1 ... +1)
% natural derivatives order:
% [d/dx, d/dy, d^2/dx^2, d^2/dy^2, d^2/dxdy]

shape = [-((eta - 1)*(xi - 1)*(eta^2 + eta + xi^2 + xi - 2))/8;
    ((eta - 1)*(xi + 1)*(eta^2 + eta + xi^2 - xi - 2))/8;
    ((eta + 1)*(xi + 1)*(- eta^2 + eta - xi^2 + xi + 2))/8;
    ((eta + 1)*(xi - 1)*(eta^2 - eta + xi^2 + xi - 2))/8;
    -((eta - 1)*(xi - 1)^2*(xi + 1))/8;
    -((eta - 1)*(xi - 1)*(xi + 1)^2)/8;
    ((eta + 1)*(xi - 1)*(xi + 1)^2)/8;
    ((eta + 1)*(xi - 1)^2*(xi + 1))/8;
    -((eta - 1)^2*(eta + 1)*(xi - 1))/8;
    ((eta - 1)^2*(eta + 1)*(xi + 1))/8;
    ((eta - 1)*(eta + 1)^2*(xi + 1))/8;
    -((eta - 1)*(eta + 1)^2*(xi - 1))/8];

naturalDerivatives(:,1) = [-((eta - 1)*(eta^2 + eta + 3*xi^2 - 3))/8;
    ((eta - 1)*(eta^2 + eta + 3*xi^2 - 3))/8;
    ((eta + 1)*(- eta^2 + eta - 3*xi^2 + 3))/8;
    -((eta + 1)*(- eta^2 + eta - 3*xi^2 + 3))/8;
    ((eta - 1)*(- 3*xi^2 + 2*xi + 1))/8;
    -((eta - 1)*(3*xi^2 + 2*xi - 1))/8;
    ((eta + 1)*(3*xi^2 + 2*xi - 1))/8;
    -((eta + 1)*(- 3*xi^2 + 2*xi + 1))/8;
    -((eta - 1)^2*(eta + 1))/8;
    ((eta - 1)^2*(eta + 1))/8;
    ((eta - 1)*(eta + 1)^2)/8;
    -((eta - 1)*(eta + 1)^2)/8];

naturalDerivatives(:,2) = [-((xi - 1)*(3*eta^2 + xi^2 + xi - 3))/8;
    -((xi + 1)*(- 3*eta^2 - xi^2 + xi + 3))/8;
    ((xi + 1)*(- 3*eta^2 - xi^2 + xi + 3))/8;
    ((xi - 1)*(3*eta^2 + xi^2 + xi - 3))/8;
    -((xi - 1)^2*(xi + 1))/8;
    -((xi - 1)*(xi + 1)^2)/8;
    ((xi - 1)*(xi + 1)^2)/8;
    ((xi - 1)^2*(xi + 1))/8;
    ((xi - 1)*(- 3*eta^2 + 2*eta + 1))/8;
    -((xi + 1)*(- 3*eta^2 + 2*eta + 1))/8;
    ((xi + 1)*(3*eta^2 + 2*eta - 1))/8;
    -((xi - 1)*(3*eta^2 + 2*eta - 1))/8];

naturalDerivatives(:,3) = [-(3*xi*(eta - 1))/4;
    (3*xi*(eta - 1))/4;
    -(3*xi*(eta + 1))/4;
    (3*xi*(eta + 1))/4;
    -((3*xi - 1)*(eta - 1))/4;
    -((3*xi + 1)*(eta - 1))/4;
    ((3*xi + 1)*(eta + 1))/4;
    ((3*xi - 1)*(eta + 1))/4;
    0;
    0;
    0;
    0];

naturalDerivatives(:,4) = [-(3*eta*(xi - 1))/4;
    (3*eta*(xi + 1))/4;
    -(3*eta*(xi + 1))/4;
    (3*eta*(xi - 1))/4;
    0;
    0;
    0;
    0;
    -((3*eta - 1)*(xi - 1))/4;
    ((3*eta - 1)*(xi + 1))/4;
    ((3*eta + 1)*(xi + 1))/4;
    -((3*eta + 1)*(xi - 1))/4];

naturalDerivatives(:,5) = [1/2 - (3*xi^2)/8 - (3*eta^2)/8;
    (3*eta^2)/8 + (3*xi^2)/8 - 1/2;
    1/2 - (3*xi^2)/8 - (3*eta^2)/8;
    (3*eta^2)/8 + (3*xi^2)/8 - 1/2;
    xi/4 - (3*xi^2)/8 + 1/8;
    1/8 - (3*xi^2)/8 - xi/4;
    (3*xi^2)/8 + xi/4 - 1/8;
    (3*xi^2)/8 - xi/4 - 1/8;
    eta/4 - (3*eta^2)/8 + 1/8;
    (3*eta^2)/8 - eta/4 - 1/8;
    (3*eta^2)/8 + eta/4 - 1/8;
    1/8 - (3*eta^2)/8 - eta/4;
    ];

end