function globalMio
  # x1' = -o*x1 + 1*g(t)
  # x2' = 1*x1 + (-o)*x2
  # x3' = 1*x2 + (-o)*x3
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
    h3(j) = 1*x2(j) + (-o)*x3(j);
    x3(j+1) = x3(j) + Dt*h3(j);

  endfor
  figure(1)
  plot(t,x3,'b')
  grid on
  title('x3(t)')

  figure(2)
  plot(tr,h3(1:N),'b')
  grid on
  title('h3(t)')

  ## 2 ## Función registro gr
  figure(3)
  plot(tr,gr,'b')
  grid on
  title('Función de registro gr(t)')

  ## 3 ## H3, TDF de h3(t)=dx3/dt. gráfico de módulo
  H3=fft(h3);

  figure(4)
  stem(abs(H3(1:70)))
  grid on
  title('H3(k), TDF de h3(t)')

  ## 4 ## Gr, TDF de gr(t)
  Gr=fft(gr);

  figure(5)
  stem(abs(Gr(1:N/2)))
  grid on
  title('Gr(k), TDF de gr(t)')

  ## 5 ## Ycf(k) = H3(k)*Gr(k)
  for k=1:N
    Ycf(k)=H3(k)*Gr(k);
  endfor

  figure(6)
  stem(abs(Ycf(1:N/20)))
  grid on
  title('Ycf(k), H3*Gr')

  ## 6 ##
  ycf=Dt*ifft(Ycf);
  figure(7)
  plot(tr,ycf,'b')
  grid on
  title('ycf(k), TDF inversa de Ycf(k)')

endfunction
