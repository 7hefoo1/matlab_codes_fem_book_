% ................................................................
% MATLAB codes for Finite Element Analysis
% problem21.m
% free vibrations of laminated plates using Q9 elements
% See reference:
% K. M. Liew, Journal of Sound and Vibration,
% Solving the vibration of thick symmetric laminates
% by Reissner/Mindlin plate theory and the p-Ritz method, Vol. 198,
% Number 3, Pages 343-360, 1996
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
h = 0.001; rho = 1; I = h^3/12;
[AMatrix,BMatrix,DMatrix,SMatrix,Q] = liewMaterial(h);

% mesh generation
L = 1;
numberElementsX = 5;
numberElementsY = 5;
numberElements=numberElementsX*numberElementsY;

[nodeCoordinates, elementNodes] = ...
    rectangularMesh(2*L,L,numberElementsX,numberElementsY,'Q9');
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

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
    formMassMatrixMindlinlaminated5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,rho,h,I,'Q9','third');

% boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC5dof('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);

% eigenproblem: free vibrations
numberOfModes = 12;
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
    stiffness,mass,numberOfModes);

omega = sqrt(eigenvalues);

% Liew, p-Ritz
D0 = Q(2,2)*h^3/12;%e2*h^3/12/(1-miu12*miu21);
% dimensionless omega
omega_bar = omega*L*L/pi/pi*sqrt(rho*h/D0);

% sort out eigenvalues
[omega,ii] = sort(omega);
modes = modes(:,ii);

% drawing mesh and deformed shape
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