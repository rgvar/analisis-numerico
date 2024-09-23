function sistemas_euler

    # y1'= -10 * y1 + 4 * y2
    # y2'= -4 * y1

    Dt= 1E-03;
    NDt=30000;
    a = 4;

    f = @ (y1,y2,t) [0.001*y1 + a*y2 + sin(2*t) * 0; -a*y1 + 0.001*y2 + sin(2*t) * 15];

    t=zeros(1,NDt);
    y=zeros(2,NDt);
    k1=zeros(2,1);

    t(1)=0;
    y(1,1)=0;
    y(2,1)=11;
    for j=1: NDt + 1
        k1(:) = Dt * f(y(1,j), y(2,j),t(j));
        t(j + 1)= t(j) + Dt;
        y(:,j + 1)= y(:,j) + k1;
        #fprintf("j=%d ; t=%f; y1=%f ; y2=%f ;k1_y1=%f ; k1_y2=%f ; \n",j,t(j),y(1,j),y(2,j),k1(1),k1(2));
    endfor

    figure(1)
    plot(y(1,:),y(2,:),'b')
    grid on

endfunction
