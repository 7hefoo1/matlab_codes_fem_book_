% ................................................................
% MATLAB codes for Finite Element Analysis
% problem21bFgm.m
% free vibrations of FGM plates using Q9 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
thickness = 17.6e-6;
n = 10; % power-law index
[AMatrix,BMatrix,DMatrix,SMatrix,Inertia] = ...s
    reddyFgmMaterial(thickness,n);

% mesh generation
L = 20*thickness;
numberElementsX = 10; numberElementsY = 10;
numberElements=numberElementsX*numberElementsY;

[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L,L,numberElementsX,numberElementsY,'Q9');
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

figure;
drawingMesh(nodeCoordinates,elementNodes,'Q9','-');
axis equal

numberNodes = size(xx,1);

% GDof: global number of degrees of freedom
GDof = 5*numberNodes;

% stiffness and mass matrices
stiffness = formStiffnessMatrixMindlinlaminated5dof ...
    (GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,AMatrix, ...
    BMatrix,DMatrix,SMatrix,'Q9','third','complete');

[mass] = ...
    formMassMatrixFgmPlate5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,Inertia,...
    'Q9','third');

% boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC5dof('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% eigenproblem: free vibrations
numberOfModes = 12;
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
    stiffness,mass,numberOfModes);

omega = sqrt(eigenvalues);

% sort out eigenvalues
[omega,ii] = sort(omega);
modes = modes(:,ii);

% dimensionless omega
omega(1)*sqrt(1.22e3*L^4/1.44e9/thickness^2)

%drawing mesh and deformed shape
modeNumber = 1;
displacements = modes(:,modeNumber);

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