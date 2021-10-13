function drawEigenmodes1D(modeNumber,numberNodes,V1,x)

for j = 1:modeNumber
    u = [V1(1:numberNodes,j)];
    subplot(modeNumber,1,j);
    plot(x',u,'k.-','markersize',12)
    set(gca,'fontsize',12)
end

end