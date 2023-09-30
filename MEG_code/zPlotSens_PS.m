function [ cz ] = zPlotSens_PS( d, dotSize, cax, cbar )
% zPlotSens
% d is the data (a vector of length 275)
% dotSize is the dot size
% cax is the color axis limits [min max], or can be [] for automatic
% cbar is a boolean for whether to include a colorbar

if nargin < 2
    dotSize = 100; cax = []; cbar = true;
    minb = min(d); maxb = max(d);
elseif nargin < 3
    cax = []; cbar = true;
    minb = min(d); maxb = max(d);
elseif nargin < 4
    cbar = true;
    if isempty(cax)
        minb = min(d); maxb = max(d);
    else
        minb = cax(1); maxb = cax(2);
    end
else
    if isempty(cax)
        minb = min(d); maxb = max(d);
    else
        minb = cax(1); maxb = cax(2);
    end
end

load('xyMEG'), 
xy(:,[71,128,132])=[];% missing MLO42 & MLT54 & MRC12 

cs = d;

figure,set(gcf,'color','white')
% if length(cs)==sum(~eyeBlinkSens)
%     scatter(xy(1,~eyeBlinkSens), xy(2,~eyeBlinkSens), dotSize, cs, 'fill')
% else
    scatter(xy(1,:), xy(2,:), dotSize, cs, 'fill')
% end

if ~isempty(cax)
    caxis(cax)
end

if cbar
    cz = colorbar;
    set(cz,'YTick',minb:((maxb-minb)/8):maxb)
    stLabels = arrayfun(@(x)num2str(x,'%0.2f'), minb:((maxb-minb)/8):maxb, 'Uniform', false);
    set(cz, 'YTickLabel', stLabels);
end

end

