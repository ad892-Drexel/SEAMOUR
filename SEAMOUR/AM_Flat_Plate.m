function [kx,ky,kz]=AM_Flat_Plate(length,width,breadth)


ar_plate = [0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8];
coeff_k = [0 0.4 0.57 0.68 0.77 0.81 0.85 0.88 0.91 0.92 0.93 0.945 0.955 0.96 0.96 0.96 0.96];
p=polyfit(ar_plate,coeff_k,8);
x=linspace(0,8,801);
rho = 997;
y=polyval(p,x);

% figure(1)
% plot(x,y)
% grid on
% xlabel('Aspect ratio of rectangular flat plate(b^2/s)')
% ylabel('Coefficient of additional mass(k)')


[xu,yv,zw] = ar(length,breadth,width);
kx=polyval(p,xu);
ky=polyval(p,yv);
kz=polyval(p,zw);

end