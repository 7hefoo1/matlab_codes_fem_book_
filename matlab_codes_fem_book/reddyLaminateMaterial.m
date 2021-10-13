function [AMatrix,BMatrix,DMatrix,SMatrix,Inertia] = ...
    reddyLaminateMaterial(thickness)
%%% REDDY TIME TRANSIENT EXAMPLE

% plate thickness
h = thickness;

stack = [0 90]; % antisymmetric cross-ply

n_lam = length(stack);
hk = h/n_lam;

% reddy orthotropic properties
E2 = 2.1e6; E1 = 25*E2;
G12 = 0.5*E2; G13 = 0.5*E2; G23 = 0.2*E2;
nu12 = 0.25; nu21 = nu12*E2/E1;
rho0 = 8e-6;
kapa = 5/6;

% Reduced stiffness constants
Q(1,1) = E1/(1-nu12*nu21);
Q(1,2) = nu12*E2/(1-nu12*nu21);
Q(2,2) = E2/(1-nu12*nu21);
Q(6,6) = G12;
Q(4,4) = G23;
Q(5,5) = G13;

A = zeros(6);
B = zeros(6);
D = zeros(6);
I0 = 0; I1 = 0; I2 = 0; % inertias
barQ = zeros(6,6,n_lam);
for k = 1:length(stack)
    theta = stack(k);
    barQ(:,:,k) = effective_props(Q,theta);
    
    zk_ = (-h/2 + hk*k);    % z_k+1
    zk = (-h/2 + hk*(k-1)); % z_k
    
    A = A +     (zk_   - zk  ).*barQ(:,:,k);
    B = B + 1/2*(zk_^2 - zk^2).*barQ(:,:,k);
    D = D + 1/3*(zk_^3 - zk^3).*barQ(:,:,k);
    
    I0 = I0 +     (zk_   - zk  ).*rho0;
    I1 = I1 + 1/2*(zk_^2 - zk^2).*rho0;
    I2 = I2 + 1/3*(zk_^3 - zk^3).*rho0;
end

AMatrix = [A(1,1),A(1,2),A(1,6);
           A(1,2),A(2,2),A(2,6);
           A(1,6),A(2,6),A(6,6)];
BMatrix = [B(1,1),B(1,2),B(1,6);
           B(1,2),B(2,2),B(2,6);
           B(1,6),B(2,6),B(6,6)];
DMatrix = [D(1,1),D(1,2),D(1,6);
           D(1,2),D(2,2),D(2,6);
           D(1,6),D(2,6),D(6,6)];
SMatrix = [kapa*A(4,4),kapa*A(4,5);
           kapa*A(4,5),kapa*A(5,5)];

Inertia = [I0 0 0 0 0; 0 I2 0 I1 0; 0 0 I2 0 I1;
    0 I1 0 I0 0; 0 0 I1 0 I0];
end

% effective properties according to the orientation theta.
function barQ = effective_props(Q,thetak)
theta = deg2rad(thetak);
cc = cos(theta);
ss = sin(theta);

barQ(1,1) = Q(1,1)*cc^4 + 2*(Q(1,2)+2*Q(6,6))*cc^2*ss^2 ...
    + Q(2,2)*ss^4;
barQ(1,2) = (Q(1,1) + Q(2,2) -4*Q(6,6))*cc^2*ss^2 ...
    + Q(1,2)*(cc^4 + ss^4);
barQ(2,2) = Q(1,1)*ss^4 + 2*(Q(1,2)+2*Q(6,6))*cc^2*ss^2 ...
    + Q(2,2)*cc^4;
barQ(1,6) = (Q(1,1) -Q(1,2) -2*Q(6,6))*cc^3*ss ...
    + (Q(1,2) - Q(2,2) +2*Q(6,6))*cc*ss^3;
barQ(2,6) = (Q(1,1) - Q(1,2) -2*Q(6,6))*cc*ss^3 ...
    + (Q(1,2) - Q(2,2) +2*Q(6,6))*cc^3*ss;
barQ(6,6) = (Q(1,1) + Q(2,2) -2*Q(1,2) -2*Q(6,6))*cc^2*ss^2 ...
    + Q(6,6)*(cc^4 + ss^4);
barQ(4,4) = Q(4,4)*cc^2 + Q(5,5)*ss^2;
barQ(4,5) = (Q(5,5) - Q(4,4))*cc*ss;
barQ(5,5) = Q(5,5)*cc^2 + Q(4,4)*ss^2;
end