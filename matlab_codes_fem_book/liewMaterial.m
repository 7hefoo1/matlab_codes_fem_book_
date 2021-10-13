function [AMatrix,BMatrix,DMatrix,SMatrix,Q] = ...
    liewMaterial(h)

kapa = pi*pi/12;
% symbolic stiffness computation
syms phi

% liew material (orthotropic)
e2=1;e1=40*e2;g23=0.5*e2;g13=0.6*e2;g12=g13;
miu12=0.25;miu21=miu12*e2/e1;factor=1-miu12*miu21;

% angles for laminate
alfas = [0,pi/2,0]; % 3 layers [0/90/0]
% upper and lower coordinates
z(1) = -(h/2); z(2) = -(h/2)+h/3; z(3) = -z(2); z(4) = -z(1);

% [Q] in 0º orientation
qbarra(1,1,1) = e1/factor;
qbarra(1,2,1) = miu21*e1/factor;
qbarra(2,1,1) = miu12*e2/factor;
qbarra(2,2,1) = e2/factor;
qbarra(3,3,1) = g12;
qbarra(4,4,1) = kapa*g23;
qbarra(5,5,1) = kapa*g13;

% transformation matrix
T = [cos(phi)^2,sin(phi)^2,-sin(2*phi),0,0; ...
    sin(phi)^2,cos(phi)^2,sin(2*phi),0,0; ...
    sin(phi)*cos(phi),-sin(phi)*cos(phi),cos(phi)^2-sin(phi)^2,0,0; ...
    0,0,0,cos(phi),sin(phi); ...
    0,0,0,-sin(phi),cos(phi)];

% [Q] reduced stiffnesses in structural axes
qBarra = T*qbarra*T.';

Qbarra = zeros(5,5,size(alfas,2));
for s = 1:size(alfas,2)
    for i = 1:5
        for j = 1:5
            Qbarra(i,j,s) = subs(qBarra(i,j,1),phi,alfas(s));
        end
    end
end
Q = Qbarra;

% equivalent stiffnesses ______________________________________________
Astiff(5,5) = 0; Bstiff(5,5) = 0; Fstiff(5,5) = 0; Istiff(5,5) = 0;
for k = 1:size(alfas,2)
    % membrane, coupling and bending
    for i = 1:3
        for j = 1:3
            Astiff(i,j) = Astiff(i,j) + Q(i,j,k)*(z(k+1)-z(k));
            Bstiff(i,j) = Bstiff(i,j) + Q(i,j,k)*(z(k+1)^2-z(k)^2)/2;
            Fstiff(i,j) = Fstiff(i,j) + Q(i,j,k)*(z(k+1)^3-z(k)^3)/3;
        end
    end
    % shear
    for i = 4:5
        for j = 4:5
            Istiff(i,j) = Istiff(i,j) + Q(i,j,k)*(z(k+1)-z(k));
        end
    end
end

% constitutive matrices
AMatrix = Astiff(1:3,1:3);
BMatrix = Bstiff(1:3,1:3);
DMatrix = Fstiff(1:3,1:3);
SMatrix = Istiff(4:5,4:5);
end