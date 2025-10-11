function customColorBar(thresh_vals, plotColorsByIx)

Ncolors = length(thresh_vals)-1;

colormap(jet(Ncolors-1));
cb = colorbar;
if plotColorsByIx
    
    
    Ncolors/(Ncolors-1)
    cb.Ticks = 1:Ncolors+1;
    seg_len = (max(thresh_vals) - min(thresh_vals))/(length(thresh_vals)-1);
else
    
end

end