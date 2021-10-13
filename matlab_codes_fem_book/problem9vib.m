% ................................................................
% MATLAB codes for Finite Element Analysis
% problem9vib.m
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% I: second moment of area
% L: length of bar
E = 1; I = 1; EI = E*I; rho = 1; A = 2.3; rhoA = rho*A;

% generation of coordinates and connectivities
numberElements = 64;
nodeCoordinates = linspace(0,1,numberElements+1)';
L = max(nodeCoordinates);
numberNodes = size(nodeCoordinates,1);
xx = nodeCoordinates(:,1);
elementNodes = zeros(numberElements,2);
for i = 1:numberElements
    elementNodes(i,1)=i; 
    elementNodes(i,2)=i+1;
end

% for structure:
%   displacements: displacement vector
%   stiffness: stiffness matrix
%   mass: mass matrix
%   GDof: global number of degrees of freedom
GDof = 2*numberNodes; 

% stiffess matrix
[stiffness,~] = ...
    formStiffnessBernoulliBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,EI,1);

% stiffess matrix
[mass] = ...
    formMassBernoulliBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,rhoA);

% boundary conditions and solution
% clamped-clamped
% fixedNodeU =[1 2*numberElements+1]'; 
% fixedNodeV =[2 2*numberElements+2]';
% simply supported-simply supported
fixedNodeU =[1 2*numberElements+1]'; fixedNodeV =[]';
% clamped at x=0
%fixedNodeU =[1]'; fixedNodeV =[2]';

prescribedDof = [fixedNodeU;fixedNodeV];

% free vibration problem
[modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
                                 stiffness,mass,0);

% natural frequencies
omega = sqrt(eigenvalues);

% exact frequencies
omega_exact(1,1) =   pi^2*sqrt(EI/rhoA/L^4);
omega_exact(2,1) = 4*pi^2*sqrt(EI/rhoA/L^4);
omega_exact(3,1) = 9*pi^2*sqrt(EI/rhoA/L^4);