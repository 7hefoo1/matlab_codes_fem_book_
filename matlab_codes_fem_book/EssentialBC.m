function [prescribedDof,activeDof] = ...
    EssentialBC(typeBC,GDof,xx,yy,nodeCoordinates,numberNodes)
% essential boundary conditions for rectangular plates
% W: transverse displamcent
% TX: rotation about y axis
% TY: rotation about x axis

switch typeBC
    case 'ssss' % simply supported plate
        fixedNodeW =find(xx==max(nodeCoordinates(:,1))|...
            xx==min(nodeCoordinates(:,1))|...
            yy==min(nodeCoordinates(:,2))|...
            yy==max(nodeCoordinates(:,2)));
        
        fixedNodeTX =find(yy==max(nodeCoordinates(:,2))|...
            yy==min(nodeCoordinates(:,2)));
        fixedNodeTY =find(xx==max(nodeCoordinates(:,1))|...
            xx==min(nodeCoordinates(:,1)));
        
    case 'cccc' % clamped plate
        fixedNodeW =find(xx==max(nodeCoordinates(:,1))|...
            xx==min(nodeCoordinates(:,1))|...
            yy==min(nodeCoordinates(:,2))|...
            yy==max(nodeCoordinates(:,2)));
        fixedNodeTX =fixedNodeW;
        fixedNodeTY =fixedNodeTX;
        
    case 'scsc'
        fixedNodeW =find(xx==max(nodeCoordinates(:,1))|...
            xx==min(nodeCoordinates(:,1))|...
            yy==min(nodeCoordinates(:,2))|...
            yy==max(nodeCoordinates(:,2)));
        
        fixedNodeTX =find(xx==max(nodeCoordinates(:,2))|...
            xx==min(nodeCoordinates(:,2)));
        fixedNodeTY=[];
        
    case 'cccf'
        fixedNodeW =find(xx==min(nodeCoordinates(:,1))|...
            yy==min(nodeCoordinates(:,2))|...
            yy==max(nodeCoordinates(:,2)));
        
        fixedNodeTX =fixedNodeW;
        fixedNodeTY =fixedNodeTX;
end

prescribedDof = [fixedNodeW; fixedNodeTX+numberNodes; ...
    fixedNodeTY+2*numberNodes];
activeDof = setdiff([1:GDof]',prescribedDof);

end