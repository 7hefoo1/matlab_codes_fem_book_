function contourField(numberElements,elementNodes,xx,yy,field,...
                      colore_max,colore_min,numeroContorni,cmax,flippedMap)

% field: to be plotted in vector form
%%
% colore_max = 2e-4;%0;
% colore_min = -4.5e-4;%-1.2e-3;
% numeroContorni = 26;
% cmax = 26;

if colore_min == -Inf
    colore_min = min(field);
end
if colore_max == Inf
    colore_max = max(field);
end

coloreField = zeros(2,2,numberElements);
XX = zeros(2,2,numberElements);
YY = zeros(2,2,numberElements);
for iel=1:numberElements
    nd=elementNodes(iel,[1 2 4 3]);
    originalField = field(nd);
    for i = 1:2
        for j = 1:2
            XX(i,j,iel) = xx(nd(i+(j-1)*2));
            YY(i,j,iel) = yy(nd(i+(j-1)*2));
            if originalField(i+(j-1)*2) >= colore_max
                coloreField(i,j,iel) = colore_max;
            elseif originalField(i+(j-1)*2) <= colore_min
                coloreField(i,j,iel) = colore_min;
            else
                coloreField(i,j,iel) = originalField(i+(j-1)*2);
            end
        end
    end
    
end

figure; hold on;
for n=1:numberElements
    [Contorno,hh] = contourf(XX(:,:,n),YY(:,:,n),coloreField(:,:,n),linspace(colore_min,colore_max,numeroContorni));
    hh.LineWidth = 1.5;
end

axis equal

if flippedMap
    colormap(flipud(jet(cmax))) % inverse jet colormap
else
    colormap(jet(cmax)) % inverse jet colormap
end
colorbar

% change color properties according to the limits of values
caxis('auto')
if colore_min == -Inf
    colore_colorbar_minimo = gca.CLim(1);
else
    colore_colorbar_minimo = colore_min;
end
if colore_max == Inf
    colore_colorbar_massimo = gca.CLim(2);
else
    colore_colorbar_massimo = colore_max;
    caxis(gca,[colore_colorbar_minimo colore_colorbar_massimo])
end
caxis(gca,[colore_colorbar_minimo colore_colorbar_massimo])

end