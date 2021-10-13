function  [mass] = ...
    formMass2Dframe(GDof,numberElements, ...
    elementNodes,numberNodes,xx,yy,rhoA)

mass = zeros(GDof); 
% computation of the system stiffness matrix
for e = 1:numberElements
  % elementDof: element degrees of freedom (Dof)
  indice = elementNodes(e,:);       
  elementDof = [indice indice+numberNodes indice+2*numberNodes];
  xa = xx(indice(2))-xx(indice(1));
  ya = yy(indice(2))-yy(indice(1));  
  length_element = sqrt(xa*xa+ya*ya);
  cosa = xa/length_element;  
  sena = ya/length_element;
  ll = length_element;
  
 L = [cosa*eye(2) sena*eye(2) zeros(2);
     -sena*eye(2) cosa*eye(2) zeros(2);
      zeros(2,4) eye(2)];
          
oneu = 1/6*[2 1;1 2];
oneu2 = 1/70*[26 9; 9 26];
oneu3 = ll/420*[22 -13; 13 -22];
oneu4 = ll^2/420*[4 -3; -3 4];

m1 = rhoA*ll*[oneu   zeros(2,4);
            zeros(2) oneu2 oneu3;
            zeros(2) oneu3' oneu4];

    mass(elementDof,elementDof) = ...
        mass(elementDof,elementDof) + L'*m1*L;
end 

