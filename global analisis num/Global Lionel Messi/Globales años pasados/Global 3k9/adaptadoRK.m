
function adaptado
  gr = load('registro_231123.txt', '-ascii');
  N=2688;
  Dt=0.02;
  t0=0;
  tf=(N-1)*Dt;
  o=1.17;
  dim=N+1;



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

% Inicialización
x(1, 1) = 0; % Inicializar x1
x(1, 2) = 0; % Inicializar x2
x(1, 3) = 0; % Inicializar x3 (nuevo estado)
t(1) = t0;     % Inicializar el tiempo
w = 0.5;       % Parámetro (peso)

% Método de Euler mejorado (Runge-Kutta de 2º Orden)
for j = 1:N
    % Calcular k1
    k1 = Dt * f_sist_general(x(j, :)', t(j), gr(j));

    % Calcular punto intermedio (Euler explícito)
    xg = x(j, :)' + k1 / (2 * w);
    tg = t(j) + Dt / (2 * w);

    % Calcular k2
    k2 = Dt * f_sist_general(xg, tg, gr(j));

    % Actualizar estado
    x(j + 1, :) = x(j, :) + (1 - w) * k1' + w * k2';
    t(j + 1) = t(j) + Dt;
endfor

  figure(2)
  plot(tr,x(1:N,3),'b')
  grid on
  title('x3(t)-EJERCICIO 2') %Ninguna grafica es correcta, es una exponencial pero al revés, esto lo comprobe con chatGpt y dice que es el resultado esperado

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

  x3_2= x(1:N,3);

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

%La respuesta de teoria es Ca es correcta y cb no correcta

endfunction



% Sistema de ecuaciones diferenciales generalizado
function [fx] = f_sist_general(z, t, gr)
  o=1.17;
    fx(1, 1) = (-o) * z(1) + gr;        % Ecuación para x1'
    fx(2, 1) = z(1) + (-o) * z(2);      % Ecuación para x2'
    fx(3, 1) = z(2) + (-o) * z(3);      % Ecuación para x3'
end
