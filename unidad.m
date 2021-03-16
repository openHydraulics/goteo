nR=25;ng=100;h0=10;k=1.25;x=0.54;CVm=0.05;leg=0.25;sg=0.5;Dr=13.6;I0r=0.05;leR=0.5;sR=1;DpR=50;I0pR=0;

%Se calcula la distribución de presión hpR en la tubería portarramal
%Se utiliza el paquete Image -> pkg load image, expresiones estadísticas para matrices.
pkg load image;
tic;

h0resum=[];
CVqresum=[];
CVqhresum=[];
qreqresum=[];
qbrutresum=[];
qnetresum=[];
qdefresum=[];

for h0=1:0.1:20

  vectorUnos=ones(1,nR);
  matrizAcum=triu(ones(nR),0);

  xpR=(sR.*vectorUnos)*matrizAcum;
  zpR=-I0pR.*xpR;

  hpR=h0.*vectorUnos;
  hpRant=0.*vectorUnos;

  %Variación manufactura
  varManuf=1+CVm.*randn(ng,nR);

  %Se calcula iterativamente la distribución de presión en el portarramal 
  while max(abs(hpR-hpRant))>1e-3;
    hpRant=hpR;
    
    %Se llama a la función ramal.m para, dado una presión en cabeza de un ramal, calcular la distribución de presiones
    for j=nR:-1:1
      %Presiones y caudales en los goteros del ramal j
      [h,q]=ramal(hpRant(j),ng,k,x,CVm,leg,sg,Dr,I0r,varManuf(:,j));
      hr(:,nR+1-j)=h;%Matriz de presiones en los goteros de la unidad
      qr(:,nR+1-j)=q;%Matriz de caudales de los goteros de la unidad
    end

    %Caudal en los tramos de la tubeŕia portarramal
    QpR=(ones(1,ng)*qr)*transpose(matrizAcum);
    hpR=h0.*vectorUnos-zpR-(0.465.*QpR.^1.75.*DpR.^-4.75.*(leR+sR))*matrizAcum;
    %++cont;
  end

  %Se muestra en pantalla el tiempo transcurrido hasta completar el cálculo completo de la distribución de presiones en la unidad para un valor de presión en cabeza de la misma
  info=[h0 toc];
  disp(info)

  qh=(k.*hr.^x);

  %Resultados
  CVq=std2(qr)/mean2(qr);
  CVqh=std2(qh)/mean2(qh);
  CVqm=sqrt(CVq^2-CVqh^2);

  h0resum=[h0resum h0];
  CVqresum=[CVqresum CVq];
  CVqhresum=[CVqhresum CVqh];
  
  qreq=mean(prctile(qr,25/2));
  qbrut=mean2(qr);
  qnet=mean2((qr>qreq).*qreq+(qr<=qreq).*qr);
  qdef=mean2((qr<qreq).*(qreq-qr));
  
  qreqresum=[qreqresum qreq];
  qbrutresum=[qbrutresum qbrut];
  qnetresum=[qnetresum qnet];
  qdefresum=[qdefresum qdef];
  
  Ra=qnetresum./qbrutresum;
  Cd=qdefresum./qreqresum;
  relHbHreq=(1-Cd)./Ra;
  
##  if h0==10
##    qr10=qr;
##    hr10=hr;
##    qreq10=mean(prctile(qr,25/2));
##    qbrut10=mean2(qr);
##    qnet10=mean2((qr10>qreq10).*qreq10+(qr10<=qreq10).*qr10);
##    qdef10=mean2((qr10<qreq10).*(qreq10-qr10));
##    Ra=qnet10/qbrut10;
##    Cd=qdef10/qreq10;
##  endif

end

CVqest=sqrt(CVqhresum.^2+CVm^2);

subplot(2,2,1)
%Se representa el coeficiente de variación del caudal CVq y el coeficiente de variación del caudal debido a la variación de presión CVqh en función de la presión en cabeza de la unidad h0
plot(h0resum,CVqresum,"+;CVq;")
grid on
hold on
xlabel('h0')
ylabel('CVq')
plot(h0resum,CVqest,"-;CVq;")
plot(h0resum,CVqhresum,"-;CVqh;")
hold off

%Se muestra en pantalla el valor de la presión en cabeza de la unidad que hace mínimo el coefciente de variación del caudal
disp("La altura de presión h0(m) que minimiza la variación del caudal debido a la presión es:")
disp(h0resum(find(CVqhresum==min(CVqhresum))))

subplot(2,2,2)
%Se representa el rendimiento de aplicación para una lámina requerida definida por la media del cuarto menor
plot(h0resum,Ra,"-;Ra;")
grid on
hold on
xlabel('h0')
ylabel('Ra, Cd')
plot(h0resum,Cd,"-;Cd;")
plot(h0resum,relHbHreq,"-;Hb/Hreq=(1-Cd)/Ra;")
hold off

subplot(2,2,3)
%Se representan los caudales neto y bruto
plot(h0resum,qreqresum,"-;qreq;")
grid on
hold on
xlabel('h0')
ylabel('qreq, qbrut')
plot(h0resum,qbrutresum,"-;qbrut;")
plot(h0resum,(qbrutresum./qreqresum),"-;qbrut/qreq;")
hold off