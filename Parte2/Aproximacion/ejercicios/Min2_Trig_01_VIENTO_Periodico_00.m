function Min2_trig_01_VIENTO_periodico_00
  clc, clear
  % g(t) es presión de VIENTO CIRSOC
  % datos
  p1=[1.0 0.9 0.7 0.35 0.0 -0.5 -0.8 -1.1 -1.23 -1.3 -1.3 -1.28 -1.25 ];
  p2=[-1.2 -1.2 -1.2 -1.2 -1.2 ];
  p3=[-1.25 -1.28  -1.3 -1.3 -1.23 -1.1 -0.8 -0.5 0.0 0.35 0.7 0.9 1.0];       
  p= [p1  p2 -1.2  p2  p3];
  
  x=0:10:360;                   %angulos en grados
  x=x*pi/180;                   %angulo en radianes
  figure (1)
   plot(x,p,'or')
   grid on
   title ('función dato')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % Inicio de Min 2  con Base {1; cos(k x) }
  N=length(x)-1 % Cantidad de Intervalos entre 0 y 2PI (puntos datos en periodo)
  m= 10          % Cantidad de Armónicos en la Base Max es 18
  g(:,1)=p(1:N);          % Datos para Min 2
  P=zeros(N,1);         % Aproximación resultante
  r=zeros(N,1);         % Residuos
  FI=zeros(N,m+1);    % Matriz FI
  alfa=zeros(m+1,1);  % Coeficientes de FIT*FI alfa =FIT*g
  
  kw=zeros(1,m+1);     % Armónicos
  Am  =zeros(1,m+1);     % Amplitud por cada Armónico
  
  %Elementos de la BASE para Min2 trigonometricos
  FI(:,1)=1;
  kw(1)=0;
  for k=1:m
     kw(k+1)= k;
     for j=1:N
      FI(j,k+1)  =cos(k*x(j));
     end
  end
  FI
  FI'*FI
  D=diag(diag(FI'*FI));
  b=FI'*g(1:N,1);
  for j=1:m+1
    alfa(j,1)=b(j)/D(j,j);
  end
  disp('los coeficientes b de MIn 2 son:')
  b 
  disp('los coeficientes alfa de MIn 2 son:')
  alfa
%
% COMPARA P de Min 2 con g(t)
% 
  P=FI*alfa;
  r=g-P;
    
  figure (2)
  plot(x(1:N), P,'-xb',x(1:N),g,'or')
  grid on
  title ('g(t); P(t) discretos')
  
  figure (3)
  plot(x(1:N), P,'-xb',x(1:N),g,'or',x(1:N),r,'*r')
  grid on
  title ('g(t); P(t); r(t) discretos')

% ARMA Rta en Frecuencia
  Am(1)=abs(alfa(1,1));
  fase(1)=0;
  for k=2:m+1
     Am(k)=abs(alfa(k,1));
    fase(k)=0;
  end
  disp('la amplitud de c para cada kw son:')
  Am
  disp('la fase para cada kw son:')
  fase
  
 
  figure (4)
  subplot(2,1,1)
  stem(kw, abs(Am),'ob')
  grid on
  title ('Amplitud para cada k-ésima frecuencia')
  %
  subplot(2,1,2)
  plot(kw, fase,'ob')
  grid on
  title ('Fase para cada k-ésima frecuencia')
  % FIN de  Rta en Frecuencia
endfunction
