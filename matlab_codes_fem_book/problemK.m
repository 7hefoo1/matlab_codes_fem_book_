% ................................................................
% MATLAB codes for Finite Element Analysis
% problemK.m
% Kirchhoff plate in bending
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% isotropic material
E = 10920; poisson = 0.30; thickness = 0.01; I = thickness^3/12;
D11 = E*thickness^3/12/(1-poisson^2);
D22 = D11; D12 = poisson*D11; D66 = (1-poisson)*D11/2;

% orthotropic material
% E1 = 31.8e6; E2 = 1.02e6;
% poisson12 = 0.31; G12 = 0.96e6;
% poisson21 = poisson12*E2/E1;
% thickness = 0.01; I = thickness^3/12;
% D11 = E1*I/(1-poisson12*poisson21);
% D22 = E2*I/(1-poisson12*poisson21);
% D12 = poisson12*D22; D66 = G12*I;


% matrix C
% bending part
C_bending = [D11 D12 0; D12 D22 0; 0 0 D66];

% load
P = -1;

% 3: non-conforming 4 node element
% 4: conforming 4 node element
dof_per_node = 3; % number of kinematic parameters per node

% mesh generation
L = 1;
numberElementsX = 20;
numberElementsY = numberElementsX;
numberElements = numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L,L,numberElementsX,numberElementsY,'Q4');
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q4','-');
axis equal

numberNodes = size(xx,1);
% GDof: global number of degrees of freedom
GDof = dof_per_node*numberNodes;

[stiffness] = formStiffnessMatrixK(GDof,numberElements, ...
            elementNodes,numberNodes,nodeCoordinates, ...
            C_bending,'complete',dof_per_node);
        
[force] = formForceVectorK(GDof,numberElements,elementNodes, ...
            numberNodes,nodeCoordinates,P,'complete',dof_per_node);

% boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% solution
displacements = solution(GDof,prescribedDof,stiffness,force);

% displacements
disp('Displacements')
jj = 1:GDof; format
f = [jj; displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)

format long
% isotropic dimensionless deflection
D1 = E*thickness^3/12/(1-poisson^2);
min(displacements(1:numberNodes))*D1/L^4

% orthotropic dimensionless deflection
% min(displacements(1:numberNodes))*(D12+2*D66)/P/L^4

% surface representation
figure; hold on
for k = 1:size(elementNodes,1)
    patch(nodeCoordinates(elementNodes(k,1:4),1),...
          nodeCoordinates(elementNodes(k,1:4),2),...
          displacements(elementNodes(k,1:4)),...
          displacements(elementNodes(k,1:4)))
end
set(gca,'fontsize',18)
view(45,45)