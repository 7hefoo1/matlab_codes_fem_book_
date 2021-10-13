function [KG] = ...
    formGeometricStiffnessMindlinlaminated5dof(GDof,...
    numberElements,elementNodes,numberNodes,...
    nodeCoordinates,sigmaMatrix,thickness,elemType,...
    quadTypeB,quadTypeS)

% computation of geometric stiffness for laminated plate element

% KG: geometric matrix
KG = zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations] = gaussQuadrature(quadTypeB);

% cycle for element
for e = 1:numberElements
    % indice: nodal connectivities for each element
    % elementDof: element degrees of freedom
    indice = elementNodes(e,:);
    elementDof = [indice indice+numberNodes indice+2*numberNodes ...
        indice+3*numberNodes indice+4*numberNodes];
    ndof = length(indice);
    
    % cycle for Gauss point
    for q = 1:size(gaussWeights,1)
        GaussPoint = gaussLocations(q,:);
        xi = GaussPoint(1);
        eta = GaussPoint(2);
        
        % shape functions and derivatives
        [shapeFunction,naturalDerivatives] = ...
            shapeFunctionsQ(xi,eta,elemType);
        
        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,invJacobian,XYderivatives]=...
            Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
        
        % geometric matrix (w)
        G_b = zeros(2,5*ndof);
        G_b(1,1:ndof)  = XYderivatives(:,1)';
        G_b(2,1:ndof)  = XYderivatives(:,2)';
        KG(elementDof,elementDof) = KG(elementDof,elementDof) + ...
            G_b'*sigmaMatrix*thickness*G_b*...
            gaussWeights(q)*det(Jacob);
        
        % geometric matrix (u)
        G_a1 = zeros(2,5*ndof);
        G_a1(1,3*ndof+1:4*ndof)  = XYderivatives(:,1)';
        G_a1(2,3*ndof+1:4*ndof)  = XYderivatives(:,2)';
        KG(elementDof,elementDof) = KG(elementDof,elementDof) + ...
            G_a1'*sigmaMatrix*thickness*G_a1*...
            gaussWeights(q)*det(Jacob);
        
        % geometric matrix (v)
        G_a2 = zeros(2,5*ndof);
        G_a2(1,4*ndof+1:5*ndof)  = XYderivatives(:,1)';
        G_a2(2,4*ndof+1:5*ndof)  = XYderivatives(:,2)';
        KG(elementDof,elementDof) = KG(elementDof,elementDof) + ...
            G_a2'*sigmaMatrix*thickness*G_a2*...
            gaussWeights(q)*det(Jacob);
        
    end  % Gauss point
    
end    % element

% shear stiffness matrix

% Gauss quadrature for shear part
[gaussWeights,gaussLocations] = gaussQuadrature(quadTypeS);

% cycle for element
for e = 1:numberElements
    % indice : nodal condofectivities for each element
    % elementDof: element degrees of freedom
    indice = elementNodes(e,:);
    elementDof = [ indice indice+numberNodes indice+2*numberNodes ...
        indice+3*numberNodes indice+4*numberNodes];
    ndof = length(indice);
    
    for q = 1:size(gaussWeights,1)
        GaussPoint = gaussLocations(q,:);
        xi = GaussPoint(1);
        eta = GaussPoint(2);
        
        % shape functions and derivatives
        [shapeFunction,naturalDerivatives] = ...
            shapeFunctionsQ(xi,eta,elemType);
        
        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,invJacobian,XYderivatives] = ...
            Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
        
        % Geometric matrix
        G_s1 = zeros(2,5*ndof);
        G_s1(1,ndof+1:2*ndof)     = XYderivatives(:,1)';
        G_s1(2,ndof+1:2*ndof)     = XYderivatives(:,2)';
        KG(elementDof,elementDof) = KG(elementDof,elementDof) + ...
            G_s1'*sigmaMatrix*thickness^3/12*G_s1*...
            gaussWeights(q)*det(Jacob);
        
        G_s2 = zeros(2,5*ndof);
        G_s2(1,2*ndof+1:3*ndof)   = XYderivatives(:,1)';
        G_s2(2,2*ndof+1:3*ndof)   = XYderivatives(:,2)';
        KG(elementDof,elementDof) = KG(elementDof,elementDof) + ...
            G_s2'*sigmaMatrix*thickness^3/12*G_s2*...
            gaussWeights(q)*det(Jacob);
        
    end  % gauss point
end    % element

end