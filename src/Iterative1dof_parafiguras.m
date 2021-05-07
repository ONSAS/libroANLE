%%%  Ejemplo Numerico 2  %%%%
%%% Solucion Iterativa %%%%
clc, clear all, close all

% Punto de Arranque
x0 = 0;

% Definicion de v
v = 1;

%Definicion de f(x)
f = @(x) x-3/2*x.^2+1/2*x.^3;

% Definicion de F_x(x)
Fx = @(x) 1-3*x+3/2*x.^2;

% Fijamos valor de Lambda_k
Lk = .19;

% Parametros de Criterios de Parada
EpsilonTol = 1e-10;    MAXITER = 20;

%%%% Solucion Iterativa con N-R %%%%
xnr(:,1)=x0;

% Comienza secuencia de iteraciones
i=1;   Error = inf;

while and(Error > EpsilonTol, i<=MAXITER)
  Deltax = -(f(xnr(:,i))-Lk*v)/Fx(xnr(:,i)) ; % Incremento N-R
  xnr(:,i+1) = xnr(:,i) + Deltax; % Paso N-R
  Error = abs(f(xnr(:,i+1))-Lk*v); % Evaluacion Error de Convergencia
  i = i+1;
end

%%%% Solucion Iterativa con N-R Modif. %%%%
xnrm(:,1)=x0;

% Comienza secuencia de iteraciones
i=1;   Error = inf;

while and(Error > EpsilonTol, i<=MAXITER)
  Deltax = -(f(xnrm(:,i))-Lk*v)/Fx(xnrm(:,1)) ; % Incremento N-R Modif
  xnrm(:,i+1) = xnrm(:,i) + Deltax; % Paso N-R Modif.
  Error = abs(f(xnrm(:,i+1))-Lk*v); % Evaluacion Error de Convergencia
  i = i+1;
end


ms=6; lw=2.2;

% Genero graficas de solucion para NR y NR Modif
subplot(1,2,1)

hold on
% Newton-Raphson
Nnr = length(xnr);
gridf = 0:.01:.5;
Solexact=f(gridf);
xnrz(1:2:2*Nnr-1) = xnr;     xnrz(2:2:2*Nnr-1) = xnr(2:Nnr);
ynrz(1:2:2*Nnr-1) = f(xnr);  ynrz(2:2:2*Nnr-1) = Lk;
plot(gridf,Solexact,'-r','markersize',ms,'linewidth',lw)
plot(xnr,f(xnr),'ok','markersize',ms,'linewidth',lw)

plot(xnrz,ynrz,'-.b','markersize',ms,'linewidth',lw*0.6)
plot([-1,6],[Lk,Lk],'k-.','markersize',ms,'linewidth',lw)
axis([0 0.5 0 0.255])

labx=xlabel('$x$'), laby=ylabel('$\lambda$')
legend('Solución Exacta', 'Solución Numérica')


set(gca, 'fontsize', 15 )

copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 12);
set(labx, "FontSize", 14);
set(laby, "FontSize", 14);



subplot(1,2,2)

hold on
% Newton-Raphson Modficado
Nnrm = length(xnrm);
gridf = 0:.01:.5;
Solexact=f(gridf);
xnrmz(1:2:2*Nnrm-1) = xnrm;     xnrmz(2:2:2*Nnrm-1) = xnrm(2:Nnrm);
ynrmz(1:2:2*Nnrm-1) = f(xnrm);  ynrmz(2:2:2*Nnrm-1) = Lk;
plot(gridf,Solexact,'-r','markersize',ms,'linewidth',lw)
plot(xnrm,f(xnrm),'ok','markersize',ms,'linewidth',lw)

plot(xnrmz,ynrmz,'-.b','markersize',ms,'linewidth',lw*0.6)
plot([-1,6],[Lk,Lk],'k-.','markersize',ms,'linewidth',lw)
axis([0 0.5 0 0.255])

labx=xlabel('$x$'), laby=ylabel('$\lambda$')
legend('Solución Exacta', 'Solución Numérica')

set(gca, 'fontsize', 15 )

copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 12);
set(labx, "FontSize", 14);
set(laby, "FontSize", 14);

print('Fig4Epscargabaja','-depslatex')
