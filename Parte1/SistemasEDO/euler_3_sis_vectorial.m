function RK2_03_con_Euler_Vectorial
% y1' =  (-10y1 + 4 y2)    con y1(0)=5
% y2' =  (-4 y1 + 0 y2)    con y2(0)=3
% yex1= (1/3)*exp(-2*t)+(14/3)*exp(-8*t));
% yex2= (2/3)*exp(-2*t)+(7/3)*exp(-8*t));
clc,clear
% Datos
y10=5;     % Valores Iniciales
y20=3;
t0=0;
Dt=0.02;  % incremento de tiempo
NDt=100;  % cantidad de Dt a realizar
% Dimensionamiento
t=zeros(1,NDt);   % vector fila para el tiempo
y=zeros(2,NDt);   % Matriz para el vector solución y(t)
k1=zeros(2,1);
k2=zeros(2,1);
% Inicialización del primer estado solución
t(1)=t0;
y(1,1)=y10;
y(2,1)=y20;
% Euler
for j=1:NDt-1
   k1=Dt*f_sist_1(y(:,j),t(j));
   y(:,j+1)  = y(:,j) + k1;
   t(j+1)    = t(j)   + Dt;
end
figure(1)
plot(t,y(1,:),'b',t,y(2,:),'r');
title ('Solución de EDO con Euler')
legend('y(1)', 'y(2)')
grid on
% Fin de EULER
%
%Comparacion con Sol Exacta
% yex1= (1/3)*exp(-2*t)+(14/3)*exp(-8*t));
% yex2= (2/3)*exp(-2*t)+(7/3)*exp(-8*t));
for i=1:NDt
    yex(i)=(1/3)*exp(-2*t(i))+(14/3)*exp(-8*t(i));
    er(i)=abs(yex(i)-y(1,i));
end
normer=norm(er,inf)

figure(3)
plot(t,er)
title ('Error en y(1) entre Sol Exacta y Solución de EDO con Euler')
grid on

figure (4)
plot (t,y(1,:),'--b',t,yex,'r' )
title ('Solución de EDO')
legend('y(1) aprox', 'y(1) exacta')
grid on
end
function [fy]=f_sist_1(z,x)
  fy(1,1)=(-10*z(1) + 4*z(2));
  fy(2,1)=(-4*z(1)+0*z(2));
end
