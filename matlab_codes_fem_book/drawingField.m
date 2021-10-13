function drawingField(coord,element,tipoElemento,field)
  
 
if ( nargin == 4 )
  nodesoff=0;
end
  
if ( size(field) == size(element) )
  elementalField=1;
else
  elementalField=0;
end

if (size(coord,2) < 3)
   for c=size(coord,2)+1:3
      coord(:,c)=[zeros(size(coord,1),1)];
   end
end

holdState=ishold;
hold on

if ( strcmp(tipoElemento,'Q9') )      % Q9 elemento
%   ord=[1,5,2,6,3,7,4,8,1];
  ord=[1,2,3,4,1];
elseif ( strcmp(tipoElemento,'Q4') )  % Q4 elemento
  ord=[1,2,3,4,1];
elseif ( strcmp(tipoElemento,'L2') )  % L2 elemento
  ord=[1,2];   
end

for e=1:size(element,1)
  
   xpt=coord(element(e,ord),1);
   ypt=coord(element(e,ord),2);      
   zpt=coord(element(e,ord),3);
   
   if ( elementalField )
     fpt=field(e,ord);
   else
     fpt=field(element(e,ord));
   end
   
   fill3(xpt,ypt,zpt,fpt)
end

shading interp
axis equal
      
if ( ~holdState )
  hold off
end
