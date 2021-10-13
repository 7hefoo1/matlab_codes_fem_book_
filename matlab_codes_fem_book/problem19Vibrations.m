% ................................................................
% MATLAB codes for Finite Element Analysis
% problem19Vibrations.m
% Mindlin plate in free vibrations Q4 elements
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear; close all

% materials
E = 10920; poisson = 0.30; G = E/2/(1+poisson); thickness = 0.1;
rho = 1; I = thickness^3/12;

% kapa = 0.8601; % cccc / cccf case
% kapa = 0.822; % scsc case
kapa = 5/6;  % ssss case

% constitutive matrix
% bending part
C_bending = I*E/(1-poisson^2)* ...
    [1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
% shear part
C_shear = kapa*thickness*E/2/(1+poisson)*eye(2);

% mesh generation
L = 1;
numberElementsX = 10;
numberElementsY = 10;
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
GDof = 3*numberNodes;

% computation of the system stiffness and mass matrices
[stiffness] = ...
    formStiffnessMatrixMindlin(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,C_shear, ...
    C_bending,'Q4','complete','reduced');

[mass]=...
    formMassMatrixMindlin(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,thickness,...
    rho,I,'Q4','complete');

% boundary conditions
[prescribedDof,activeDof] = ...
    EssentialBC('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% free vibrations
numberModes = 15;
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
    stiffness,mass,numberModes);

omega = sqrt(eigenvalues);

% sort out eigenvalues
[omega,ii] = sort(omega);
modes = modes(:,ii);

% dimensionless omega
omega_bar = omega*L*sqrt(rho/G)

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

% plot of the first 4 eigenmodes
for k = 1:4
    modeNumber = k;
    displacements = modes(:,modeNumber);
    contourField(numberElements,elementNodes,xx,yy,...
        displacements(1:numberNodes),Inf,-Inf,11,22,1)
    title(['\omega(',num2str(k),') = ',num2str(omega(k))])
    colorbar off; box on;
    set(gca,'linewidth',5,'fontsize',18)
end