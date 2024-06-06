function [cmap] = get_cmap(model)

blue = [0 0 255]./256;
green = [0.4660 0.6740 0.1880];
white = [200 200 200]./256; %grey
red = [255 0 0]./256;
cyan = [0 255 255]./256;
magenta = [0.4940 0.1840 0.5560];
black = [0 0 0];
yellow = [255 255 0]./256;
dk_blue = [0 0.4470 0.7410];
orange = [0.9290 0.6940 0.1250];
cherry = [0.6350 0.0780 0.1840];

if model == 2
    cmap = {blue, cyan, green, dk_blue, orange,red, cherry, magenta};
else
    cmap = {blue, cherry, green, orange, red, cyan, magenta, dk_blue};
end
    