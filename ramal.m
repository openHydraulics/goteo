%Cálculo de la distribución de caudal en un ramal de riego por goteo
h0=10;n=100;k=1.25;x=0.54;CVm=0.05;le=0.5;s=0.5;D=13.6;I0=0.05;

vectorUnos=ones(1,n);
matrizAcum=tril(ones(n),0);

xR=(s.*vectorUnos)*transpose(matrizAcum);
zR=-I0.*xR;

%Variación manufactura
varManuf=1+CVm.*randn(1,n);

h=h0.*vectorUnos-zR;
q=(k.*h.^x).*varManuf;

%Resultados
h0
CVq=std(q)/mean(q)
CVqh=std(k.*h.^x)/mean(k.*h.^x)
CVqm=sqrt(CVq^2-CVqh^2)

plot(xR,zR,"-;z(x);");hold on;plot(xR,h,"-;h(x);");plot(xR,-q,"o;q(x);");hold off;