% ................................................................
function  [mass] = ...
    formMassBernoulliBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx,rhoA)

mass = zeros(GDof);
% calculation of the system mass matrix
for e = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [2*(indice(1)-1)+1 2*(indice(2)-1) ...
        2*(indice(2)-1)+1 2*(indice(2)-1)+2];
    % length of element
    LElem = xx(indice(2))-xx(indice(1));
    k1 = rhoA*LElem/420*[156   22*LElem 54 -13*LElem;
                        22*LElem 4*LElem^2 13*LElem -3*LElem^2;
                        54 13*LElem 156 -22*LElem;
                        -13*LElem -3*LElem^2 -22*LElem 4*LElem^2];
    
    % mass matrix
    mass(elementDof,elementDof) = ...
        mass(elementDof,elementDof) + k1;
end

end
