function [modes,eigenvalues] = eigenvalue(GDof,prescribedDof,...
    stiffness,mass,maxEigenvalues)
% function to find solution in terms of global displacements
% GDof: number of degree of freedom
% prescribedDof: bounded boundary dofs
% stiffness: stiffness matrix
% mass: mass matrix
% maxEigenvalues: maximum eigenvalues to be computed. If 0 all the
% eigenvalues are requested (suggested for beam structures)

%%
activeDof = setdiff((1:GDof)', prescribedDof);
if maxEigenvalues == 0
    [V,D] = eig(stiffness(activeDof,activeDof),...
        mass(activeDof,activeDof));
else
    [V,D] = eigs(stiffness(activeDof,activeDof),...
        mass(activeDof,activeDof),maxEigenvalues,'smallestabs');
end

eigenvalues = diag(D);
modes = zeros(GDof,length(eigenvalues));
modes(activeDof,:) = V;

end