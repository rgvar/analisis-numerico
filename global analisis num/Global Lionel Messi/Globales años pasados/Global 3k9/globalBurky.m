function globalBurky
  gr = load('registro_231123.txt', '-ascii');
  N=2688;
  Dt=0.02;
  t0=0;
  tf=(N-1)*Dt;
  o=1.17;
  dim=N+1;

  ## 1 ## Euler explícito graficar x3 y h3
  x1(1)=0;
  x2(1)=0;
  x3(1)=0;
  t(1)=0;

  for j=1:N
    tr(j)=t0+(j-1)*Dt;
  endfor

  figure(1)
  plot(tr,gr(1:N),'b')
  grid on
  title('Registro-EJERCICIO 1') %Figura 11 correcta

  for j=1:N
    % función escalón unitaria
    if tr(j)>=0
      g(j)=1;
    else
      g(j)=0;
    endif
  endfor

  for j=1:N # Euler
    t(j+1) = t(j) + Dt;
    x1(j+1) = x1(j) + Dt*((-o)*x1(j) + 1*g(j));
    x2(j+1) = x2(j) + Dt*(1*x1(j) + (-o)*x2(j));
    x3(j+1) = x3(j) + Dt*(1*x2(j) + (-o)*x3(j));
  endfor

  figure(2)
  plot(tr,x3(1:N),'b')
  grid on
  title('x3(t)-EJERCICIO 2')

for i=1:N
    h(i)=exp(-o*tr(i));
  endfor

  %Gráfico de la función h(t)
  figure(3)
  plot(tr,h(1:N),'-b')
  grid on
  title('Función h(t)-EJERCICIO 3') % Figura 31 correcta
  legend('Función h(t)')

  %------Convolución-------

  %h3(t) asociada a x3(t)
% ej 4
   h3_1=Dt*conv(h,h);
  h3=Dt*conv(h,h3_1);

  figure(4)
  plot(tr,h3(1:N),'b')
  grid on
  title('h3(t)-EJERCICIO 4') % Figura 42 correcta


  X3 = Dt*conv(h3,gr);
  figure(5)
  plot(tr,X3(1:N),'b')
  grid on
  title('x3(t) entre h3 y registro - Ejercicio 5') % Figura 52 correcta


  %Ejercicio 6

  x3_2= x3;

  x3_5=X3;

  %Integral con Trapecios Múltiple
  I1=0;
  for k=1:N-1
    I1= I1 + (x3_2(k)^2 + x3_2(k+1)^2)*(Dt/2); %48.769
  endfor
  disp('El valor de la Integral I1 es: ')
  I1

  Id=0;
  for k=1:N-1
    Id= Id + ((x3_2(k)-x3_5(k))^2 + (x3_2(k+1)-x3_5(k+1))^2)*(Dt/2);  %0.3309
  endfor
  disp('El valor de la Integral Id es: ')
  Id
  printf("El valor de la Integral Id es: %.6f\n", Id);

%La respuesta de teoria es Ca es correcta y cb no correcta

endfunction
