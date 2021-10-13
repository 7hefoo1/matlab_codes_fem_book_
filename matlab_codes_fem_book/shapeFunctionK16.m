function [shape,naturalDerivatives] = shapeFunctionK16(xi,eta)
% shape function and derivatives for conforming elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi and eta
% xi, eta: natural coordinates (-1 ... +1)
% natural derivatives order:
% [d/dx, d/dy, d^2/dx^2, d^2/dy^2, d^2/dxdy]

shape = [((eta - 1)^2*(eta + 2)*(xi - 1)^2*(xi + 2))/16;
    -((eta - 1)^2*(eta + 2)*(xi + 1)^2*(xi - 2))/16;
    ((eta + 1)^2*(eta - 2)*(xi + 1)^2*(xi - 2))/16;
    -((eta + 1)^2*(eta - 2)*(xi - 1)^2*(xi + 2))/16;
    ((eta - 1)^2*(eta + 2)*(xi - 1)^2*(xi + 1))/16;
    ((eta - 1)^2*(eta + 2)*(xi - 1)*(xi + 1)^2)/16;
    -((eta + 1)^2*(eta - 2)*(xi - 1)*(xi + 1)^2)/16;
    -((eta + 1)^2*(eta - 2)*(xi - 1)^2*(xi + 1))/16;
    ((eta - 1)^2*(eta + 1)*(xi - 1)^2*(xi + 2))/16;
    -((eta - 1)^2*(eta + 1)*(xi + 1)^2*(xi - 2))/16;
    -((eta - 1)*(eta + 1)^2*(xi + 1)^2*(xi - 2))/16;
    ((eta - 1)*(eta + 1)^2*(xi - 1)^2*(xi + 2))/16;
    ((eta - 1)^2*(eta + 1)*(xi - 1)^2*(xi + 1))/16;
    ((eta - 1)^2*(eta + 1)*(xi - 1)*(xi + 1)^2)/16;
    ((eta - 1)*(eta + 1)^2*(xi - 1)*(xi + 1)^2)/16;
    ((eta - 1)*(eta + 1)^2*(xi - 1)^2*(xi + 1))/16];

naturalDerivatives(:,1) = [(3*(xi^2 - 1)*(eta - 1)^2*(eta + 2))/16;
    -(3*(xi^2 - 1)*(eta - 1)^2*(eta + 2))/16;
    (3*(xi^2 - 1)*(eta + 1)^2*(eta - 2))/16;
    -(3*(xi^2 - 1)*(eta + 1)^2*(eta - 2))/16;
    -((eta - 1)^2*(eta + 2)*(- 3*xi^2 + 2*xi + 1))/16;
    ((eta - 1)^2*(eta + 2)*(3*xi^2 + 2*xi - 1))/16;
    -((eta + 1)^2*(eta - 2)*(3*xi^2 + 2*xi - 1))/16;
    ((eta + 1)^2*(eta - 2)*(- 3*xi^2 + 2*xi + 1))/16;
    (3*(xi^2 - 1)*(eta - 1)^2*(eta + 1))/16;
    -(3*(xi^2 - 1)*(eta - 1)^2*(eta + 1))/16;
    -(3*(xi^2 - 1)*(eta - 1)*(eta + 1)^2)/16;
    (3*(xi^2 - 1)*(eta - 1)*(eta + 1)^2)/16;
    -((eta - 1)^2*(eta + 1)*(- 3*xi^2 + 2*xi + 1))/16;
    ((eta - 1)^2*(eta + 1)*(3*xi^2 + 2*xi - 1))/16;
    ((eta - 1)*(eta + 1)^2*(3*xi^2 + 2*xi - 1))/16;
    -((eta - 1)*(eta + 1)^2*(- 3*xi^2 + 2*xi + 1))/16];

naturalDerivatives(:,2) = [(3*(eta^2 - 1)*(xi - 1)^2*(xi + 2))/16;
    -(3*(eta^2 - 1)*(xi + 1)^2*(xi - 2))/16;
    (3*(eta^2 - 1)*(xi + 1)^2*(xi - 2))/16;
    -(3*(eta^2 - 1)*(xi - 1)^2*(xi + 2))/16;
    (3*(eta^2 - 1)*(xi - 1)^2*(xi + 1))/16;
    (3*(eta^2 - 1)*(xi - 1)*(xi + 1)^2)/16;
    -(3*(eta^2 - 1)*(xi - 1)*(xi + 1)^2)/16;
    -(3*(eta^2 - 1)*(xi - 1)^2*(xi + 1))/16;
    -((xi - 1)^2*(xi + 2)*(- 3*eta^2 + 2*eta + 1))/16;
    ((xi + 1)^2*(xi - 2)*(- 3*eta^2 + 2*eta + 1))/16;
    -((xi + 1)^2*(xi - 2)*(3*eta^2 + 2*eta - 1))/16;
    ((xi - 1)^2*(xi + 2)*(3*eta^2 + 2*eta - 1))/16;
    -((xi - 1)^2*(xi + 1)*(- 3*eta^2 + 2*eta + 1))/16;
    -((xi - 1)*(xi + 1)^2*(- 3*eta^2 + 2*eta + 1))/16;
    ((xi - 1)*(xi + 1)^2*(3*eta^2 + 2*eta - 1))/16;
    ((xi - 1)^2*(xi + 1)*(3*eta^2 + 2*eta - 1))/16];

naturalDerivatives(:,3) = [(3*xi*(eta - 1)^2*(eta + 2))/8;
    -(3*xi*(eta - 1)^2*(eta + 2))/8;
    (3*xi*(eta + 1)^2*(eta - 2))/8;
    -(3*xi*(eta + 1)^2*(eta - 2))/8;
    ((3*xi - 1)*(eta - 1)^2*(eta + 2))/8;
    ((3*xi + 1)*(eta - 1)^2*(eta + 2))/8;
    -((3*xi + 1)*(eta + 1)^2*(eta - 2))/8;
    -((3*xi - 1)*(eta + 1)^2*(eta - 2))/8;
    (3*xi*(eta - 1)^2*(eta + 1))/8;
    -(3*xi*(eta - 1)^2*(eta + 1))/8;
    -(3*xi*(eta - 1)*(eta + 1)^2)/8;
    (3*xi*(eta - 1)*(eta + 1)^2)/8;
    ((3*xi - 1)*(eta - 1)^2*(eta + 1))/8;
    ((3*xi + 1)*(eta - 1)^2*(eta + 1))/8;
    ((3*xi + 1)*(eta - 1)*(eta + 1)^2)/8;
    ((3*xi - 1)*(eta - 1)*(eta + 1)^2)/8];

naturalDerivatives(:,4) = [(3*eta*(xi - 1)^2*(xi + 2))/8;
    -(3*eta*(xi + 1)^2*(xi - 2))/8;
    (3*eta*(xi + 1)^2*(xi - 2))/8;
    -(3*eta*(xi - 1)^2*(xi + 2))/8;
    (3*eta*(xi - 1)^2*(xi + 1))/8;
    (3*eta*(xi - 1)*(xi + 1)^2)/8;
    -(3*eta*(xi - 1)*(xi + 1)^2)/8;
    -(3*eta*(xi - 1)^2*(xi + 1))/8;
    ((3*eta - 1)*(xi - 1)^2*(xi + 2))/8;
    -((3*eta - 1)*(xi + 1)^2*(xi - 2))/8;
    -((3*eta + 1)*(xi + 1)^2*(xi - 2))/8;
    ((3*eta + 1)*(xi - 1)^2*(xi + 2))/8;
    ((3*eta - 1)*(xi - 1)^2*(xi + 1))/8;
    ((3*eta - 1)*(xi - 1)*(xi + 1)^2)/8;
    ((3*eta + 1)*(xi - 1)*(xi + 1)^2)/8;
    ((3*eta + 1)*(xi - 1)^2*(xi + 1))/8];

naturalDerivatives(:,5) = [(9*(eta^2 - 1)*(xi^2 - 1))/16;
    -(9*(eta^2 - 1)*(xi^2 - 1))/16;
    (9*(eta^2 - 1)*(xi^2 - 1))/16;
    -(9*(eta^2 - 1)*(xi^2 - 1))/16;
    -(3*(eta^2 - 1)*(- 3*xi^2 + 2*xi + 1))/16;
    (3*(eta^2 - 1)*(3*xi^2 + 2*xi - 1))/16;
    -(3*(eta^2 - 1)*(3*xi^2 + 2*xi - 1))/16;
    (3*(eta^2 - 1)*(- 3*xi^2 + 2*xi + 1))/16;
    -(3*(xi^2 - 1)*(- 3*eta^2 + 2*eta + 1))/16;
    (3*(xi^2 - 1)*(- 3*eta^2 + 2*eta + 1))/16;
    -(3*(xi^2 - 1)*(3*eta^2 + 2*eta - 1))/16;
    (3*(xi^2 - 1)*(3*eta^2 + 2*eta - 1))/16;
    ((- 3*eta^2 + 2*eta + 1)*(- 3*xi^2 + 2*xi + 1))/16;
    -((- 3*eta^2 + 2*eta + 1)*(3*xi^2 + 2*xi - 1))/16;
    ((3*eta^2 + 2*eta - 1)*(3*xi^2 + 2*xi - 1))/16;
    -((3*eta^2 + 2*eta - 1)*(- 3*xi^2 + 2*xi + 1))/16];

end