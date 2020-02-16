nR=25;ng=100;h0=10;k=1.25;x=0.54;CVm=0.05;leg=0.75;sg=0.5;Dr=13.6;I0r=0.05;leR=1.5;sR=1;DpR=50;I0pR=0;

%Se calcula la distribución de presión hpR en la tubería portarramal
%se utiliza el paquete Image -> pkg load image
pkg load image;
tic;

h0resum=[];
CVqresum=[];
CVqhresum=[];

for h0=1:0.1:20

  vectorUnos=ones(1,nR);
  matrizAcum=triu(ones(nR),0);

  xpR=(sR.*vectorUnos)*matrizAcum;
  zpR=-I0pR.*xpR;

  hpR=h0.*vectorUnos;
  hpRant=0.*vectorUnos;

  %Variación manufactura
  varManuf=1+CVm.*randn(ng,nR);

  cont=0;

  while max(abs(hpR-hpRant))>1e-3;
    hpRant=hpR;
    for j=nR:-1:1
      %Presiones y caudales en los goteros del ramal j
      [h,q]=ramal(hpRant(j),ng,k,x,CVm,leg,sg,Dr,I0r,varManuf(:,j));
      hr(:,nR+1-j)=h;%Matriz de presiones en los goteros de la unidad
      qr(:,nR+1-j)=q;%Matriz de caudales de los goteros de la unidad
    end

    %Caudal en los tramos de la tubeŕia portarramal
    QpR=(ones(1,ng)*qr)*transpose(matrizAcum);
    hpR=h0.*vectorUnos-zpR-(0.465.*QpR.^1.75.*DpR.^-4.75.*(leR+sR))*matrizAcum;
    ++cont;
  end

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

end

CVqest=sqrt(CVqhresum.^2+CVm^2);

plot(h0resum,CVqresum,"+;CVq;")
grid on
hold on
xlabel('h0')
ylabel('CVq')
plot(h0resum,CVqest,"-;CVq;")
plot(h0resum,CVqhresum,"-;CVqh;")

disp("La altura de presión h0(m) que minimiza la variación del caudal debido a la presión es:")
disp(h0resum(find(CVqhresum==min(CVqhresum))))

%surf(hr)
%hold on
%surf(qr)