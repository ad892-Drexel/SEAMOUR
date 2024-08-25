function [xu,yv,zw]=ar(length,thickness,width)
    xu = length/thickness;
    yv = width/thickness;
    zw = length/width;
end