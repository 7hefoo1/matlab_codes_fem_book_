function [force] = ...
    formForceVectorMindlin5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,P,...
    elemType,quadType)

% computation of force vector for laminated plate element

% force: force vector
force = zeros(GDof,1);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations] = gaussQuadrature(quadType);

% cycle for element
for e = 1:numberElements
    % indice: nodal connectivities for each element
    indice = elementNodes(e,:);
    
    % cycle for Gauss point
    for q = 1:size(gaussWeights,1)
        GaussPoint = gaussLocations(q,:);
        GaussWeight = gaussWeights(q);
        xi = GaussPoint(1);
        eta = GaussPoint(2);
        
        % shape functions and derivatives
        [shapeFunction,naturalDerivatives] = ...
                    shapeFunctionsQ(xi,eta,elemType);
        
        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,invJacobian,XYderivatives] = ...
            Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
        
        % force vector
        force(indice) = force(indice) + ...
            shapeFunction*P*det(Jacob)*GaussWeight;
    end  % end Gauss point loop
    
end  % end element loop

end