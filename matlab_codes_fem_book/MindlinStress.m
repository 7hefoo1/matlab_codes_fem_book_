function [stress,shear] = MindlinStress(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,displacements,...
    C_shear,C_bending,h,elemType,quadTypeB,quadTypeS)
%% Mindlin Stress
% computes normal and shear stresses according to Mindlin
% theory note that transverse shear stresses are not corrected


%% normal stresses
% 1: sigma_xx
% 2: sigma_yy
% 3: tau_xy
stress = zeros(numberElements,4,3);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations] = gaussQuadrature(quadTypeB);

% cycle for element
for e = 1:numberElements
    % indice : nodal connectivities for each element
    % indiceB: element degrees of freedom
    indice = elementNodes(e,:);
    indiceB = [indice indice+numberNodes indice+2*numberNodes];
    nn = length(indice);
    
    % cycle for Gauss point
    for q = 1:size(gaussWeights,1)
        pt = gaussLocations(q,:);
        wt = gaussWeights(q);
        xi = pt(1);
        eta = pt(2);
        
        % shape functions and derivatives
        [shapeFunction,naturalDerivatives] = ...
                    shapeFunctionsQ(xi,eta,elemType);
        
        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,invJacobian,XYderivatives] = ...
            Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
        
        % [B] matrix bending
        B_b = zeros(3,3*nn);
        B_b(1,nn+1:2*nn)        = XYderivatives(:,1)';
        B_b(2,2*nn+1:3*nn)      = XYderivatives(:,2)';
        B_b(3,nn+1:2*nn)        = XYderivatives(:,2)';
        B_b(3,2*nn+1:3*nn)      = XYderivatives(:,1)';
        
        % stresses
        strain = h/2*B_b*displacements(indiceB);
        stress(e,q,:) = C_bending*strain;
        
    end  % Gauss point
end    % element


%% shear stresses
% 1: tau_xz
% 2: tau_yz
% by constitutive equations
shear = zeros(numberElements,1,2);

% Gauss quadrature for shear part
[gaussWeights,gaussLocations] = gaussQuadrature(quadTypeS);

% cycle for element
for e = 1:numberElements
    % indice : nodal connectivities for each element
    % indiceB: element degrees of freedom
    indice = elementNodes(e,:);
    indiceB = [ indice indice+numberNodes indice+2*numberNodes];
    nn = length(indice);
    
    % cycle for Gauss point
    for q = 1:size(gaussWeights,1)
        pt = gaussLocations(q,:);
        wt = gaussWeights(q);
        xi = pt(1);
        eta = pt(2);
        
        % shape functions and derivatives
        [shapeFunction,naturalDerivatives] = ...
                    shapeFunctionsQ(xi,eta,elemType);
        
        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,invJacobian,XYderivatives] = ...
            Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
        
        % [B] matrix shear
        B_s = zeros(2,3*nn);
        B_s(1,1:nn)       = XYderivatives(:,1)';
        B_s(2,1:nn)       = XYderivatives(:,2)';
        B_s(1,nn+1:2*nn)  = shapeFunction;
        B_s(2,2*nn+1:3*nn)= shapeFunction;
        
        sliding = B_s*displacements(indiceB);
        shear(e,q,:) = C_shear*sliding;
        
    end  % end gauss point loop
end    % end element loop


end