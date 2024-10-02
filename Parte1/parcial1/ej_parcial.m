function ej_parcial
clc,clear
% Datos
y10=0;
y20=0;
t0=0;
fs=5000
Dt=1/fs;
NDt=fs*4;
w=1/2; % Euler mejorado
% Inicialización variables
t=zeros(1,NDt);
y=zeros(2,NDt);
b=load("p1_3k9_01.txt","-ascii");
k1=zeros(2,1);
k2=zeros(2,1);
% Inicialización del primer estado
t(1)=t0;
y(1,1)=y10;
y(2,1)=y20;
% Runge-Kutta 2do orden
for j=1:NDt-1
   k1=Dt*f_sist_1(y(:,j),t(j),b(j));
     tg=t(j)+Dt/(2*w);
     yg=y(:,j)+k1/(2*w);
   k2= Dt*f_sist_1(yg,tg,b(j));
   y(:,j+1)= y(:,j) + (1-w)*k1 + w*k2;
   t(j+1)= t(j) + Dt;
end
figure(1)
plot(t,y(1,:),'r');
title ('Solución de EDO con Runge-Kutta')
grid on

end
function [fy]=f_sist_1(z,x,b)
  fy(1,1)=(0*z(1) + 1*z(2));
  fy(2,1)=(-14400*z(1)+(-24)*z(2)) + b;
end
