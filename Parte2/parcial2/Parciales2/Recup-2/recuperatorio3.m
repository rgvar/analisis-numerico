function recuperatorio3

%%===============DATOS=====================
gn=load('registro-13-nov.txt','-ascii');
N= length(gn)
dt=0.083776;
tn=zeros(1,N);

for j=1:N
  tn(j)=dt*(j-1);
end
figure (1)
plot (tn,gn,'r')
grid on
legend("Funcion g(tn)")
%==========NORMA CUADRATICA=================
Norma_2_g=0;
for k=1:N
Norma_2_g = Norma_2_g + gn(k)^2;
end
Norma_2_g
%===========TDF====================

Tp=dt*N;
dw=(2*pi)/Tp

G_tdf= fft(gn,N);
for k=1:N
  mod_G(k)=abs(G_tdf(k));
end

Norma_2_G=0;
for k=1:N
Norma_2_G= Norma_2_G + mod_G(k) * mod_G(k);
end
Norma_2_G

 %Crear un vector de frecuencias
frecuencias = dw * (0:(N-1));

  % Graficar el módulo de la TDF
  figure(2)
  #grafico reducido
  stem(abs(G_tdf(1:N/9)))
  title("Transformada de Fourier (modulo)")
%=============CONVOLUCION============
A1=1/2;
kc=8;
a=0.3;
p=kc*dw;
for i=1:N
    h(i)= A1 *tn(i) * exp(-p*tn(i));
end

yc=dt*conv(h,gn);
figure(3)
plot(tn,yc(1:N),'-b')
grid on
legend("Convolucion")

figure(4)
plot(tn,h(1:N),'b')
grid on
legend("h")

figure (5)
plot(tn,gn,'r',tn,yc(1:N),'b')
grid on
legend ("G(t), Yc")
end

