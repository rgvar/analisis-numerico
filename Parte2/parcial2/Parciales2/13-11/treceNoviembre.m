function treceNoviembre
  # función g con datos brindados
  gn = load("registro-13-nov.txt", "-ascii");
  # datos iniciales
  Dt=0.083776;
  t0=0;
  N=length(gn) # para saber cantidad de datos, son 1024
 # vector tn (abscisas)
  for k=1:N
    tn(k)=t0+(k-1)*Dt;
  endfor
%Ejercicio 1
  Norma_2_g=norm(gn,2) %El 2 porque es raiz cuadrada
    figure(1)
  plot(tn,gn,'r')
  title('Generación de la función discreta dato')
  grid on
  title('figura 1')
  %Ejercicio 2
  Tp = Dt*N;
  Dw= (2*pi)/Tp
  G_tdf=fft(gn,N);
  for k=1:N
      mod_G(k)=abs(G_tdf(k));
  endfor
  Norma_2_G = norm(mod_G)
  figure(2);
  stem(mod_G(1:N/9),'b');
  title('Grafica del modulo');
    grid on;
  %Ejercicio3
  kc=8;
  A1=1/3;
  a=0.3;
  w=kc*Dw;
  p=a*w;
 for i=1:N
    h(i)= A1 * exp(-p*tn(i)) * sin(w * tn(i))/w;
end
  yc=Dt*conv(h,gn);
  figure(3)
  plot(tn,yc(1:N),'b',tn,gn(1:N),'-r')
  grid on;
  legend('conv','g(t)')

  figure(4)
  plot(tn,h(1:N),'b')
    grid on;
endfunction
