colormap default
    cmap = colormap(hot);
    cmap = cmap([1:44,64],:);
    cmap = flipud(cmap);
    colormap(cmap)