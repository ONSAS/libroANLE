%% Ejemplo Numerico 1 : Solucion Incremental %%
clc, clear all, close all

% Condiciones Iniciales
x0 = 0; L0 = 0;

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

x(:,1) = x0 ;   L(1)=0;

% Defino Incremento Delta Lambda
DL = .02;

% Comienza secuencia de incrementos
i=1;
while L(end)<Lmax 
  x(:,i+1) = x(:,i) + DL * h(x(:,i)); % Forward Euler
  L(i+1) = L(i)+DL;
  i = i+1;
end

% Genero graficas de las soluciones
subplot(1,2,1)  % Solucion Por debajo de punto limite
plot(x(1,1:end-2),L(1:end-2),'*-k',0:.01:.5,f(0:.01:.5),'-r')
xlabel('x'), ylabel('$\lambda$')
legend('Solucion Numerica','Solucion Exacta')

subplot(1,2,2)  % Solucion Por encima de punto limite
plot(x(1,:),L,'*-k',0:.01:2.25,f(0:.01:2.25),'-r')
xlabel('x'), ylabel('$\lambda$')
legend('Solucion Numerica','Solucion Exacta')
