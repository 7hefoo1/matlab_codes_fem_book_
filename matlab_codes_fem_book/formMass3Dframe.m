function  [mass] = ...
    formMass3Dframe(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,rho,A,Iz,Iy)

mass = zeros(GDof);
% computation of the system mass matrix
for e = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [6*indice(1)-5 6*indice(1)-4 6*indice(1)-3 ...
        6*indice(1)-2 6*indice(1)-1 6*indice(1)...
        6*indice(2)-5 6*indice(2)-4 6*indice(2)-3 ...
        6*indice(2)-2 6*indice(2)-1 6*indice(2)] ;
    x1 = nodeCoordinates(indice(1),1);
    y1 = nodeCoordinates(indice(1),2);
    z1 = nodeCoordinates(indice(1),3);
    x2 = nodeCoordinates(indice(2),1);
    y2 = nodeCoordinates(indice(2),2);
    z2 = nodeCoordinates(indice(2),3);
    
    L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + ...
        (z2-z1)*(z2-z1));
    
    p = (Iz+Iy)/A;
    
    % lumped mass matrix
%     m = rho*A*L/2*diag([1 1 1 p 0 0 1 1 1 p 0 0]);
    
    % consistent mass matrix
    m = rho*A*L/420*[140 0 0 0 0 0 70 0 0 0 0 0;
                     0  156 0 0 0 22*L 0 54 0 0 0 -13*L;
                     0 0 156 0 -22*L 0 0 0 54 0 13*L 0;
                     0 0 0 140*p 0 0 0 0 0 70*p 0 0;
                     0 0 -22*L 0 4*L^2 0 0 0 -13*L 0 -3*L^2 0;
                     0 22*L 0 0 0 4*L^2 0 13*L 0 0 0 -3*L^2;
                     70 0 0 0 0 0 140 0 0 0 0 0;
                     0 54 0 0 0 13*L 0 156 0 0 0 -22*L;
                     0 0 54 0 -13*L 0 0 0 156 0 22*L 0;
                     0 0 0 70*p 0 0 0 0 0 140*p 0 0;
                     0 0 13*L 0 -3*L^2 0 0 0 22*L 0 4*L^2 0;
                     0 -13*L 0 0 0 -3*L^2 0 -22*L 0 0 0 4*L^2];
    
    if x1 == x2 && y1 == y2
        if z2 > z1
            Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
        else
            Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
        end
    else
        CXx = (x2-x1)/L;
        CYx = (y2-y1)/L;
        CZx = (z2-z1)/L;
        D = sqrt(CXx*CXx + CYx*CYx);
        CXy = -CYx/D;
        CYy = CXx/D;
        CZy = 0;
        CXz = -CXx*CZx/D;
        CYz = -CYx*CZx/D;
        CZz = D;
        Lambda = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];
        
    end
    R = [Lambda zeros(3,9); zeros(3) Lambda zeros(3,6);
        zeros(3,6) Lambda zeros(3);zeros(3,9) Lambda];
    
    mass(elementDof,elementDof) = ...
        mass(elementDof,elementDof) + R'*m*R;
end

end
