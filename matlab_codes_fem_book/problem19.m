% ................................................................
% MATLAB codes for Finite Element Analysis
% problem19.m
% Mindlin plate in bending Q4 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
E = 10920; poisson = 0.30; kapa = 5/6; thickness = 0.1; 
I = thickness^3/12;

% constitutive matrix
% bending part
C_bending = I*E/(1-poisson^2)* ...
    [1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
% shear part
C_shear = kapa*thickness*E/2/(1+poisson)*eye(2);

% load
P = -1;

% mesh generation
L = 1;
numberElementsX = 20;
numberElementsY = 20;
numberElements = numberElementsX*numberElementsY;
%
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L,L,numberElementsX,numberElementsY,'Q4');
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q4','-');
axis equal

numberNodes = size(xx,1);
% GDof: global number of degrees of freedom
GDof = 3*numberNodes;

% computation of the system stiffness matrix and force vector
[stiffness] = ...
    formStiffnessMatrixMindlin(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,C_shear, ...
    C_bending,'Q4','complete','reduced');

[force] = ...
    formForceVectorMindlin(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,P,'Q4','reduced');

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
D1 = E*thickness^3/12/(1-poisson^2);
min(displacements(1:numberNodes))*D1/L^4

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

% post-computation
[stress,shear] = MindlinStress(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,displacements,...
    C_shear,C_bending,thickness,'Q4','complete','reduced');