function repaso_global
  clc, clear;
  %Cargo la función g de un .txt
  g = load("registroglobal.txt", "-ascii");
  N = length(g)

  %N=2688;
  %Datos
  Dt=0.02;
  omega=1.17;
  t0=0;
  tf=N*Dt
  t=zeros(1,N);

  for i=1:N
    t(i)=t0 + Dt*(i-1);
  endfor

  figure(1)
  plot(t, g, 'b')
  grid on
  legend("g(t)")

  %------Método de Euler explícito-------
  % Dimensionamiento
  x=zeros(3,N);   % Matriz para el vector solucion x(t)

  xa=zeros(3,1);
  k1=zeros(3,1);
  % Inicializacion del primer estado solucion
  x(1,1)=0;
  x(2,1)=0;
  x(3,1)=0;

  % Euler explícito
  for j=1:N-1
    xa=x(:,j);
    ta=t(j);

    k1(1,1)= Dt* (-omega*xa(1) + 0.0*xa(2) + 0.0*xa(3) + 1*g(j));
    k1(2,1)= Dt* (1.0*xa(1) - omega*xa(2) + 0.0*xa(3) + 0.0*g(j));
    k1(3,1)= Dt* (0.0*xa(1) + 1.0*xa(2) - omega*xa(3) + 0.0*g(j));

    x(:,j+1)=xa + k1;
    t(j+1)=ta + Dt;
  end

  %Gráfico de la función x3(t)
  figure(2)
  plot(t,x(3,:),'-b')
  grid on
  legend('x3(t)')

  %------funcion h(t)------

  for i=1:N
    h(i)=exp(-omega*t(i));
  endfor

  %Gráfico de la función h(t)
  figure(3)
  plot(t,h,'-b')
  grid on
  title('Función h(t)')
  legend('Función h(t)')

  %------Convolución-------

  %h3(t) asociada a x3(t)
  h3_1=Dt*conv(h,h);
  h3=Dt*conv(h,h3_1);

  %Gráfico de la función h3(t)
  figure(4)
  plot(t,h3(1:N),'-b')
  grid on
  title('Función h3(t)')
  legend('Función h3(t)')

##  X3=Dt*conv(h3,x(2,:));
  X3=Dt*conv(h3,g);

  %Gráfico de la función x3(t)
  figure(5)
  plot(t,X3(1:N),'-b')
  grid on
  legend('x3(t)')

  %Ejercicio 6
  x3_2= x(3,:);
  x3_5=X3;

  %Integral con Trapecios Múltiple
  I1=0;
  for k=1:N
    I1= I1 + (x3_2(k)^2 + x3_2(k)^2)*(Dt/2); %48.766
  endfor
  disp('El valor de la Integral I1 es: ')
  I1

  Id=0;
  for k=1:N
    Id= Id + ((x3_2(k)-x3_5(k))^2 + (x3_2(k)-x3_5(k))^2)*(Dt/2);
  endfor
  disp('El valor de la Integral Id es: ')
  Id

endfunction
