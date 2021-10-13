% ................................................................
% MATLAB codes for Finite Element Analysis
% problem9buk.m
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% I: second moment of area
% L: length of bar
E = 1; I = 1; EI=E*I;

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
%   stability: geometric matrix
%   GDof: global number of degrees of freedom
GDof = 2*numberNodes; 

% stiffess matrix
[stiffness,~] = ...
    formStiffnessBernoulliBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx,EI,1);

% stability matrix
[stability] = ...
    formStabilityBernoulliBeam(GDof,numberElements, ...
    elementNodes,numberNodes,xx);

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
                                 stiffness,stability,0);

% natural frequencies
N0 = eigenvalues;

% exact frequencies simply-supported beam
N0_exact(1,1) =   pi^2*EI/L^2;
N0_exact(2,1) = 4*pi^2*EI/L^2;
N0_exact(3,1) = 9*pi^2*EI/L^2;