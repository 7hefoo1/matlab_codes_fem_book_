% ................................................................
% MATLAB codes for Finite Element Analysis
% problem9.m
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% I: second moment of area
% L: length of bar
E = 1; I = 1; EI=E*I;

% generation of coordinates and connectivities
numberElements = 2;
nodeCoordinates = linspace(0,1,numberElements+1)';
L = max(nodeCoordinates);
numberNodes = size(nodeCoordinates,1);
xx = nodeCoordinates(:,1);
elementNodes = zeros(numberElements,2);
for i = 1:numberElements
    elementNodes(i,1)=i; 
    elementNodes(i,2)=i+1;
end

% distributed load
P = -1;

% for structure:
%   displacements: displacement vector
%   force : force vector
%   stiffness: stiffness matrix
%   GDof: global number of degrees of freedom
GDof = 2*numberNodes; 

% stiffess matrix and force vector
[stiffness,force] = ...
    formStiffnessBernoulliBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,EI,P);

% boundary conditions and solution
% clamped-clamped
% fixedNodeU =[1 2*numberElements+1]'; 
% fixedNodeV =[2 2*numberElements+2]';
% simply supported-simply supported
fixedNodeU =[1 2*numberElements+1]'; fixedNodeV =[]';
% clamped at x=0
%fixedNodeU =[1]'; fixedNodeV =[2]';

prescribedDof = [fixedNodeU;fixedNodeV];

% solution
displacements = solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    GDof,prescribedDof)

% reordering displacements and rotations
W = displacements(1:2:2*numberNodes);
R = displacements(2:2:2*numberNodes);

% drawing nodal displacements
figure
plot(nodeCoordinates,W,'ok','markersize',8,'linewidth',1.5)
set(gca,'fontsize',18)

% graphical representation with interpolation for each element
drawInterpolatedBeam