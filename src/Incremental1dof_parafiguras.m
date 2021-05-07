%%% Ejemplo Numerico 1   %%%%
%%% Solucion Incremental %%%%
clc
clear

% Condiciones Iniciales
x0 = 0;
L0 = 0;

% Definicion de v
v = 1;

% Definicion de f(x)
f = @(x) x-3/2*x.^2+1/2*x.^3;

% Definicion de F_x(x)
Fx = @(x) 1-3*x+3/2*x.^2;

% Definicion de h(x,lambda)
h = @(x) Fx(x)\v;

%%%% Solucion Incremental con Forward Euler %%%%
% Defino hasta que Lambda incrementamos
Lmax = .25;

x(:,1)=x0;
L(1)=0;

% Defino Incremento Delta Lambda
DL = .02;

% Comienza secuencia de incrementos
i=1;
while L(end)<Lmax 
  x(:,i+1) = x(:,i) + DL * h(x(:,i)); % Forward Euler
  L(i+1) = L(i)+DL;
  i = i+1;
end

% Se grafica
ms=6; lw=2.2;

% Genero graficas de las soluciones
subplot(1,2,1)

hold on
% Solucion Por debajo de punto limite
plot(x(1,1:end-4),L(1:end-4),'o-k','markersize',ms,'linewidth',lw)
plot(0:.01:.5,f(0:.01:.5),'-r','markersize',ms,'linewidth',lw)
%~ title('Ejemplo solucion incremental')
labx=xlabel('$x$')
laby=ylabel('$\lambda$')
legend("Solución Num\\'erica",'Solución Exacta','location','SouthEast')



set(gca, 'fontsize', 15 )

copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 12);
set(labx, "FontSize", 14);
set(laby, "FontSize", 14);



%~ axisLim = axis(); axisLim(4) = 0.25;
%~ axis=axisLim


subplot(1,2,2)
hold on

% Solucion Por encima de punto limite
plot(x(1,:),L,'o-k','markersize',ms,'linewidth',lw)
plot(0:.01:2.25,f(0:.01:2.25),'-r','markersize',ms,'linewidth',lw)

%~ title('Problemas para superar punto limite')
labx=xlabel('$x$')
laby=ylabel('$\lambda$')
legend('Solución Numérica','Solución Exacta','location','NorthWest')


set(gca, 'fontsize', 15 )

copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 12);
set(labx, "FontSize", 14);
set(laby, "FontSize", 14);

print('Fig1Epscargabaja','-depslatex')
