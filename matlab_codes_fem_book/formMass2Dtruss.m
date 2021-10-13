function [mass] = ...
    formMass2Dtruss(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,xx,yy,rhoA)

mass=zeros(GDof);

% computation of the system stiffness matrix
for e = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice(1)*2-1 indice(1)*2 ...
                  indice(2)*2-1 indice(2)*2] ;
    xa = xx(indice(2))-xx(indice(1));
    ya = yy(indice(2))-yy(indice(1));
    length_element = sqrt(xa*xa+ya*ya);
    % consistent mass matrix
    k1 = rhoA*length_element/6* ...
        [2 0 1 0; 0 2 0 1;
        1 0 2 0; 0 1 0 2];
    % lumped mass matrix
%     k1 = rhoA*length_element/2*eye(4);
    mass(elementDof,elementDof) = ...
        mass(elementDof,elementDof)+k1;
end

end