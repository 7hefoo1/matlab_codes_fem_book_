% ................................................................
% MATLAB codes for Finite Element Analysis
% problem18a.m
% 2D problem: beam in bending using Q8 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
E = 10e7; poisson = 0.30;

% matrix C
C = E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

% load
P = 1e6;

% mesh generation
Lx = 5;
Ly = 1;
numberElementsX = 20;
numberElementsY = 10;
numberElements = numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(Lx,Ly,numberElementsX,numberElementsY,'Q8');
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q8','-');
axis equal

numberNodes = size(xx,1);
% GDof: global number of degrees of freedom
GDof = 2*numberNodes; 

% calculation of the system stiffness matrix
stiffness = formStiffnessMass2D(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,C,1,1,'Q8','complete');

% boundary conditions 
fixedNodeX = find(nodeCoordinates(:,1)==0);  % fixed in XX
fixedNodeY = find(nodeCoordinates(:,1)==0);  % fixed in YY
prescribedDof = [fixedNodeX; fixedNodeY+numberNodes];

% force vector (distributed load applied at xx=Lx)
force = zeros(GDof,1);
rightBord = find(nodeCoordinates(:,1)==Lx);
force(rightBord(1:2:end)+numberNodes) = P*Ly/numberElementsY/3;
force(rightBord(2:2:end)+numberNodes) = P*Ly/numberElementsY*2/3;
force(rightBord(1)+numberNodes) = P*Ly/numberElementsY/6;
force(rightBord(end)+numberNodes) = P*Ly/numberElementsY/6;

% solution
displacements = solution(GDof,prescribedDof,stiffness,force);

% displacements and deformed shape
disp('Displacements')
jj = 1:GDof; format
f = [jj; displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)
UX = displacements(1:numberNodes);
UY = displacements(numberNodes+1:GDof);
scaleFactor = 0.1;

% deformed shape
figure
drawingField(nodeCoordinates+scaleFactor*[UX UY], ...
    elementNodes,'Q9',UX);%U XX
hold on
drawingMesh(nodeCoordinates+scaleFactor*[UX UY], ...
    elementNodes,'Q8','-');
drawingMesh(nodeCoordinates,elementNodes,'Q8','--');
colorbar
title('Displacement field u_x (on deformed shape)')
axis off

% stresses at nodes
[stress,strain] = stresses2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,displacements,...
    C,'Q8','complete');

% drawing stress fields on deformed shape
figure; hold on;
drawingField2(nodeCoordinates,elementNodes,...
    scaleFactor,UX,UY,'Q4',stress(:,:,1))
axis equal
drawingMesh(nodeCoordinates,elementNodes,'Q8','--');
colorbar
title('Stress field \sigma_{xx} (on deformed shape)')
axis off

% stress extrapolation
stressExtr = zeros(numberElements,8,3);
for e = 1:numberElements
    for i = 1:3
        stressExtr(e,:,i) = [1+sqrt(3)/2 -0.5 1-sqrt(3)/2 -0.5;
           -0.5 1+sqrt(3)/2 -0.5 1-sqrt(3)/2;
           1-sqrt(3)/2 -0.5 1+sqrt(3)/2 -0.5;
           -0.5 1-sqrt(3)/2 -0.5 1+sqrt(3)/2;
           (1+sqrt(3))/4 (1+sqrt(3))/4 (1-sqrt(3))/4 (1-sqrt(3))/4;
           (1-sqrt(3))/4 (1+sqrt(3))/4 (1+sqrt(3))/4 (1-sqrt(3))/4;
           (1-sqrt(3))/4 (1-sqrt(3))/4 (1+sqrt(3))/4 (1+sqrt(3))/4;
           (1+sqrt(3))/4 (1-sqrt(3))/4 (1-sqrt(3))/4 (1+sqrt(3))/4]*...
            [stress(e,1,i);stress(e,2,i);stress(e,3,i);stress(e,4,i)];
    end
end

% stress averaging at nodes
stressAvg = zeros(numberNodes,3);
for i = 1:3
    currentStress = stressExtr(:,:,i);
    for n = 1:numberNodes
        idx = find(n==elementNodes);
        stressAvg(n,i) = sum(currentStress(idx))/length(currentStress(idx));
    end
end

% surface representation
figure; hold on
for k = 1:size(elementNodes,1)
    patch(nodeCoordinates(elementNodes(k,1:4),1),...
          nodeCoordinates(elementNodes(k,1:4),2),...
          nodeCoordinates(elementNodes(k,1:4),1)*0,...
          stressAvg(elementNodes(k,1:4),1))
end
axis equal; axis off
colorbar
title('Averaged nodal stress field \sigma_{xx}')