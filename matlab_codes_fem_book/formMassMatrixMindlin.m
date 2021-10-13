function [mass] = ...
    formMassMatrixMindlin(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,thickness,...
    rho,I,elemType,quadType)

% computation of  mass matrix for Mindlin plate element

% mass: mass matrix
mass = zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations] = gaussQuadrature(quadType);

% cycle for element
for e = 1:numberElements
    % indice : nodal connectivities for each element
    indice = elementNodes(e,:);
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
        [Jacob,invJacobian,XYderivatives] = ...
            Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
        
        % mass matrix
        mass(indice,indice) = mass(indice,indice) + ...
            shapeFunction*shapeFunction'*thickness* ...
            rho*gaussWeights(q)*det(Jacob);
        mass(indice+numberNodes,indice+numberNodes) = ...
            mass(indice+numberNodes,indice+numberNodes) + ...
            shapeFunction*shapeFunction'*I* ...
            rho*gaussWeights(q)*det(Jacob);
        mass(indice+2*numberNodes,indice+2*numberNodes) = ...
            mass(indice+2*numberNodes,indice+2*numberNodes) + ...
            shapeFunction*shapeFunction'*I* ...
            rho*gaussWeights(q)*det(Jacob);
        
    end % end Gauss point loop
end % end element loop

end