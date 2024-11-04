function Aprox_por_minimos_cuadrados
  clc, clear
  %Vector g
  g = load("parcial_2_3k10_2023.txt", "-ascii");
  global N = length(g)
  global w0=(2*pi)/N;
  global Dt = 0.12566; %Dt = TP/N
  global m =fix((N-1)/2);
  global dim=1+2*m;
  A1 =1/3;

  %Vector t
  t=zeros(N,1);
  for i=1 : N
      t(i) = (i-1) * Dt;
  endfor

  figure(1) %Grafico de gt
  plot(t,g)
  grid on

  %Entrega de valores maximos
  [gm,i]=max(g(1:N/3))
  tm = t(i,1)

  %TDF
  g_tdf = fft(g,N);
  figure(2)
  bar (abs(g_tdf))
  grid on

  g_tdf_abs = abs(g_tdf);
  [G_tdf_MAX,i]=max(g_tdf_abs)

  Valor = 0;
  kc = 0;
  for i=1:10
    if g_tdf_abs(i) > 0.05 * G_tdf_MAX
      if i > kc
        Valor = g_tdf_abs(i);
        kc = i;
      endif
    endif
  endfor
  display(kc)
  display(Valor)

  %Aprox_por_minimos_cuadrados

  phi = calcular_matriz_phi(N, dim); % Llamada a función y devuelve phi

  % PHI' * PHI * alfa = PHI' * g
  A = phi' * phi;  % MATRIZ DIAGONAL
  B = phi' * g;    % Vector de terminos independientes
  % A alfa = B



  %Vector de incognitas (vector alfa)
  alfa= zeros(dim,1); %Matriz de incognitas  [a0,a1,b1,a2,b2,...]
  % a0= B(1,1)/N, a1 = B(2,1)/(N/2),  b1 = B(3,1)/(N/2),.....
  alfa(1,1)= B(1,1)/N;
  for j=2:dim
    alfa(j,1)=B(j,1)/(N/2);
  endfor


  %Calcular P(t)
  N2 = N;
  dt2 = Dt;
  dw = w0/ Dt;
  for k=1:N2
    t2(1,k)=dt2*(k-1);
  endfor
  P=zeros(N2,1);

  for i=1: N2
    valor = alfa(1,1);
    for k=1:kc
      valor = valor + alfa(k*2,1)*cos(k*dw*(t2(1,i)-t2(1,1))) + alfa(k*2+1,1)*sin(k*dw*(t2(1,i)-t2(1,1)));
    endfor
    P(i,1)= valor;
  endfor

  figure(3)
  plot(t,P)
  grid on

 [P_Max,i]=max(P)

 %Convolucion
   p = kc * dw

   h=zeros(N,1);
  for i = 1:N
    h(i,1) = A1 *exp(-p * t(i,1));
  endfor

  y = Dt * conv(h,g);  %con 2N-1 de longitud


 [y_Max,i]=max(y)

  figure(4)
  plot(t, g, '-r', t, P, '-b',t,y(1:N),'-g')
  grid on
endfunction

function [phi] = calcular_matriz_phi(N, dim)

  phi = zeros(N, dim);
  global m;
  global w0;

  for n=1:N %primero recorro por filas
    phi(n,1)=1;
    for k=1:m %recorro por columnas pero no por dim sino por m
      phi(n,2*k)=cos(k*w0*(n-1));
      phi(n,2*k+1)=sin(k*w0*(n-1));
    endfor
    endfor
end



