%%  Ejemplo Numerico 2 - Solucion Iterativa %%
clc, clear all, close all

% Punto de Arranque
x0 = 0;

% Definicion de v
v = 1;

%Definicion de f(x) y F_x(x)
f  = @(x) x-3/2*x.^2+1/2*x.^3 ;
Fx = @(x) 1-3*x+3/2*x.^2      ;

% Fijamos valor de Lambda_k
Lk = .19;

% Parametros de Criterios de Parada
EpsilonTol = 1e-10;    MAXITER = 20;

%%%% Solucion Iterativa con N-R %%%%
xnr(:,1) = x0;    i=1;   Error = inf;
while and(Error > EpsilonTol, i<=MAXITER) % Comienzan iteraciones
  Deltax = -(f(xnr(:,i))-Lk*v)/Fx(xnr(:,i)) ; % Incremento N-R
  xnr(:,i+1) = xnr(:,i) + Deltax; % Paso N-R
  Error = abs(f(xnr(:,i+1))-Lk*v); % Evaluacion Error de Convergencia
  i = i+1;
end

%%%% Solucion Iterativa con N-R Modif. %%%%
xnrm(:,1) = x0 ;    i=1;   Error = inf;
while and(Error > EpsilonTol, i<=MAXITER) % Comienzan iteraciones
  Deltax = -(f(xnrm(:,i))-Lk*v)/Fx(xnrm(:,1)) ; % Incremento N-R Modif
  xnrm(:,i+1) = xnrm(:,i) + Deltax; % Paso N-R Modif.
  Error = abs(f(xnrm(:,i+1))-Lk*v); % Evaluacion Error de Convergencia
  i = i+1;
end

% Genero graficas de solucion para NR y NR Modif
subplot(1,2,1)   % Newton-Raphson
Nnr  = length(xnr);
grid = 0:.01:.5;
Solexact = f(grid);
xnrz(1:2:2*Nnr-1) = xnr;     xnrz(2:2:2*Nnr-1) = xnr(2:Nnr);
ynrz(1:2:2*Nnr-1) = f(xnr);  ynrz(2:2:2*Nnr-1) = Lk;
plot(grid,Solexact,'-r',xnrz,ynrz,'-.b',xnr,f(xnr),'ok',[-1,6],[Lk,Lk],'k-.')
axis([0 0.5 0 0.255])
xlabel('x'), ylabel('\lambda'), legend('Solucion Exacta')

subplot(1,2,2)  % Newton-Raphson Modficado
Nnrm = length(xnrm);
grid = 0:.01:.5;
Solexact = f(grid);
xnrmz(1:2:2*Nnrm-1) = xnrm;     xnrmz(2:2:2*Nnrm-1) = xnrm(2:Nnrm);
ynrmz(1:2:2*Nnrm-1) = f(xnrm);  ynrmz(2:2:2*Nnrm-1) = Lk;
plot(grid,Solexact,'-r',xnrmz,ynrmz,'-.b',xnrm,f(xnrm),'ok',[-1,6],[Lk,Lk],'k-.')
axis([0 0.5 0 0.255])
xlabel('x'), ylabel('\lambda'), legend('Solucion Exacta')
