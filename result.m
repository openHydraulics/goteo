%Representaciones de diversos resultados de riego y de ajustes a la distribución Normal
%Debe ejecutarse después de unidad.m

pkg load statistics;

resRiego=[];

f=[0.01:0.01:0.99];

for CVq=0:0.01:0.2
  Hi=(1-CVq.*norminv(f,0,1));
  Hr=mean(Hi(floor(3/4*numel(Hi)):numel(Hi)));
  
##  if CVq==0.07 | CVq==0.1
##    plot(f,Hi)
##    axis ij
##    axis([0 1 0 1.2])
##    grid on
##    xlabel('f')
##    ylabel('Hi/Hmed')
##    hold on
##    plot(f,Hr)
##  end
  
  Hp=mean((Hi>Hr).*(Hi-Hr));
  Hd=mean((Hi<Hr).*(Hr-Hi));
  Hn=Hr-Hd;
  Hb=Hn+Hp;
  Ra=Hn/Hb;
  Cd=Hd/Hr;
  resRiego=[resRiego; CVq Ra Cd (1-Cd)/Ra Hr];
end

##plot(resRiego(:,1),resRiego(:,2),"-;Ra;")
##grid on
##xlabel('CV')
##ylabel('Ra, Cd, (1-Cd)/Ra')
##hold on
##plot(resRiego(:,1),resRiego(:,3),"-;Cd;")
##plot(resRiego(:,1),resRiego(:,4),"-;(1-Cd)/Ra;")

%Representación distribución de caudales de emisores a partir de la matriz de caudal procedente de unidad.m y del ajuste a la distribución Normal
f2=[0:1/(ng*nR):1-1/(ng*nR)];
A=sort(reshape(qr10,1,(ng*nR)),"descend");
B=norminv(sort(f2,"descend"),mean(A),std(A));
plot(f2,A,'o')
axis([0 1 0 5])
axis ij
grid on
hold on
plot(f2,B)
xlabel('f')
ylabel('q(L/h)')

##surf(hr10)
##hold on
##surf(qr10)