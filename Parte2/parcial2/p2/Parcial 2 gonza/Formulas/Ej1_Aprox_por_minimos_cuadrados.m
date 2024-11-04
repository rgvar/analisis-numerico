function Aprox_por_minimos_cuadrados
  clc, clear
  Tp = 10;
  Td = 0;
  global N = 32;
  global w0=(2*pi)/N;
  global Dt = Tp/N;
  global m =fix((N-1)/2);
  global dim=1+2*m;

  %Vector t
  %t = linspace(0, Tp, N);
  t=zeros(N,1);
  for i=1 : N
      t(i) = (i-1) * Dt;
  endfor

  %Vector g
  g=zeros(N,1);
  for k=1:N
    if k <= N/2
      g(k,1)=1;
    else
      g(k,1)=0;
    endif
  endfor

  figure(1) %Grafico de gt
  plot(t,g, "marker", "o", "markerEdgeColor", "b", ...
  "markersize", 4, "color","red")
  grid on

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
  disp(alfa)

  %Amplitud
  [k, Amp] = calcular_amplitudes(alfa); % Llamada a función
  figure (2)  %Grafico de Amplitud con frecuencia
  bar(k,Amp,'ob')
  grid on
  legend('Amplitud para cada k-ésima frecuencia')

  %Calcular P(t)
  N2 = 100;
  Tp2 = 20;
  dt2 = Tp2/N2;
  dw = w0/ Dt;
  for k=1:N2
    t2(1,k)=dt2*(k-1);
  endfor
  P=zeros(N2,1);

  for i=1: N2
    valor = alfa(1,1);
    for k=1:m
      valor = valor + alfa(k*2,1)*cos(k*dw*(t2(1,i)-t2(1,1))) + alfa(k*2+1,1)*sin(k*dw*(t2(1,i)-t2(1,1)));
    endfor
    P(i,1)= valor;
  endfor

  figure(3)
  plot(t,g,'xb',t2,P, "marker", "o", "markerEdgeColor", "b", ...
  "markersize", 4, "color","red")
  grid on
  display(A)
  display(B)
  display(alfa)

  %Vector R(t)
  r = calcular_matriz_resto(g, P); % Llamada a función y devuelve r
  figure(4)
  plot(t,r, "marker", "o", "markerEdgeColor", "b", ...
  "markersize", 4, "color","red")
  grid on

  %norm(r,2)


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

function [k, Amp] = calcular_amplitudes(alfa) % Devuelve el vector k para poder graficar y el vector Amp que es el vector de amplitudes
    global m;

    k=zeros(1,m+1);%%vector K de frecuencias
    for j=1:m+1
      k(1,j)=j-1;
    endfor

    Amp=zeros(m+1, 1);  % vector de Amplitudes [Ao,A1,A2,...]
    Amp(1,1) = alfa(1,1); % Amp(0) = a0
    for j=1:m
      Amp(j+1,1) = sqrt((alfa(2*j,1))^2 + (alfa(2*j+1,1))^2);
    endfor
end
function [r] = calcular_matriz_resto(g, P)
  global N;
  r=zeros(N,1);
  for k=1: N
    r(k) = g(k) - P(k);
  endfor
end

