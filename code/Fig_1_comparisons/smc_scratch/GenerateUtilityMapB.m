function [X,Y,G] = GenerateUtilityMapB(opt,addnoise)
xdel=1;%resolution in x
ydel=1;%resolution in y
xRange=opt.DomainBounds.xmin:xdel:opt.DomainBounds.xmax-xdel;
yRange=opt.DomainBounds.ymin:ydel:opt.DomainBounds.ymax-ydel;

[X,Y] = meshgrid(xRange,yRange);

if (addnoise==1)
    n=0.05;
else
    n=0;
end

G = ones(22500,1);
G=max(G,0); %crop below 0
G=G./max(G); %normalize
G = G./sum(sum(G));

end

