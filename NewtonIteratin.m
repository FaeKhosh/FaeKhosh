clear all;  clc;  close all;
clear variables; clear all;  clc;  close all;
h=0.5;
g=1.5;
ai=0.355028;
aip=-0.258819;
R=100;
d=2;
tol = 0.00001; % Convergence tolerance!
K0=((3*(d-h))./(4*g*(h^3))).^(1/3); % initial value 
varr=0.01; 
tolerance = 0.000001; %7 digit accuracy is desired
maxIterations = 100; %Don't allow the iterations to continue indefinitely
for n = 1 : maxIterations
    K=K0;
    y = i.*(d-h) - ((4.*(K.^3).*(h^3).*g)./3) - ( (3.*(K.^2).*(h.^3).*aip+(4*(h.^2).*((i.*K).^(5/3)).*ai)).*((i^(1/3).*(K.^(4/3))+(3*aip)))  +  ((K.^3).*(h.^3).*aip+((3/2).*(i^(5/3)).*(h^2).*(K.^(8/3)).*ai)).*((4/3).*((i.*K).^(1/3)))  )./((i.^(1/3)).*(K.^(4/3))+(3*aip)).^2 ; 
    K=K0+ varr ;
    fp =  i.*(d-h) - ((4.*(K.^3).*(h^3).*g)./3) - ( (3.*(K.^2).*(h.^3).*aip+(4*(h.^2).*((i.*K).^(5/3)).*ai)).*((i^(1/3).*(K.^(4/3))+(3*aip)))  +  ((K.^3).*(h.^3).*aip+((3/2).*(i^(5/3)).*(h^2).*(K.^(8/3)).*ai)).*((4/3).*((i.*K).^(1/3)))  )./((i.^(1/3)).*(K.^(4/3))+(3*aip)).^2 ; 
    yprime= (fp - y)./varr;
    K1 = K0 - (y./yprime); %Do Newton's computation
    if(abs(K1 - K0) <= tolerance) %If the result is within the desired tolerance
        break; %Done, so leave the loop
    end
    K0 = K1 %Update K0 to start the process again
    x=real(K0);
    y=imag(K0);
    XX=linspace(x-0.5,x+0.5,100);
    YY=linspace(y-0.5,y+0.5,100);
    [X,Y] = meshgrid(XX,YY);
    a=-(i.*(X+(i.*Y)).*h)-(g.*((X+(i.*Y)).^4).*(h^3))./3 -(((X+(i.*Y)).^3).*(h.^3).*aip +((3/2).*(i^(5/3)).*(h^2).*((X+(i.*Y)).^(8/3)).*ai) )./((i^(1/3)).*((X+(i.*Y)).^(4/3))+(3*aip));
    Nk= a + i.*(X+ i*Y).* d;
    z=real(Nk);
    z=imag(Nk);
end
contour(XX,YY,z,150)