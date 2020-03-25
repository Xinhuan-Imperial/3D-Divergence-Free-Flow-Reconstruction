function fun=colored_quiver(x,y,z,u,v,w,scale,tit)
% draw 2D or 3D colored vector figures,and the color is determined by
% magnitude of vector, using the current colormap

fun=quiver3(x,y,z,u,v,w,scale);
set(fun, 'MaxHeadSize',1) ;
mags = sqrt(sum(cat(2,fun.UData(:),fun.VData(:),fun.WData(:)).^2, 2));
% currentColormap = colormap(winter(512));
currentColormap = colormap(parula(512));
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
set(fun.Head,'ColorBinding', 'interpolated','ColorData', reshape(cmap(1:3,:,:), [], 4).');
set(fun.Tail,'ColorBinding', 'interpolated','ColorData', reshape(cmap(1:2,:,:), [], 4).');
title(tit);
daspect([1 1 1]);
xlabel('x');ylabel('y');zlabel('z');
c=colorbar;
tmp=linspace(min(mags),max(mags),11)';

% tmp=round(tmp);
% tmp=linspace(min(mags),9,11)';
% tmp=linspace(min(mags),28,11)';
% tmp=linspace(0,20,11)';
title(c,'velocity [mm/s]');
c.TickLabels=[];
for i=1:length(tmp)
    c.TickLabels=[c.TickLabels;cellstr(num2str(tmp(i)))];
end

% view([1 1 0]);
end