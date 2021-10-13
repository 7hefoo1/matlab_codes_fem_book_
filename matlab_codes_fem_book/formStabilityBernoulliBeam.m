% ................................................................
function  [stability] = ...
    formStabilityBernoulliBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx)

stability = zeros(GDof);
% calculation of the system mass matrix
for e = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [2*(indice(1)-1)+1 2*(indice(2)-1) ...
        2*(indice(2)-1)+1 2*(indice(2)-1)+2];
    % length of element
    LElem = xx(indice(2))-xx(indice(1));
    k1 = 1/30/LElem*[36   3*LElem -36 3*LElem;
                     3*LElem 4*LElem^2 -3*LElem -LElem^2;
                     -36 -3*LElem 36 -3*LElem;
                     3*LElem -LElem^2 -3*LElem 4*LElem^2];
    
    % stability matrix
    stability(elementDof,elementDof) = ...
        stability(elementDof,elementDof) + k1;
end

end
