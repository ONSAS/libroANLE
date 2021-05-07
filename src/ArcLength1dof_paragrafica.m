%%%  Ejemplo Numerico 3  %%%%
%%% Solucion con Metodo de Longitud de Arco %%%%
clc, clear all, close all

% Punto Inicial
x0 = [0];     L0 = 0;

% Definicion de Ecuacion No-Lineal
v = 1;
f = @(x) x-3/2*x.^2+1/2*x.^3;
F_x = @(x) 1-3*x+3/2*x.^2;

% Pasos de Solucion con ArcLength
Nsteps = 50;

x = x0;   L = L0;

% Parametros de la restriccion de Long. Arco
psi = 1;    Dl = .083;

for k=1:Nsteps 
   % Inicializacion de Incrementos para Iteracion. F.Euler, 
   % Se ajusta el signo de DL en base a pasaje de un punto critico.
   DL_ki = Dl/2;
   
   if F_x(x(k))>0 
      DL_ki = DL_ki;      % si la rama es ascendente usa +Forw.Euler
   else
      DL_ki = -DL_ki;     % si la rama es ascendente usa -Forw.Euler
   end  
   Dx_ki = DL_ki*v/F_x(x(k)); % Incremento Forw.Euler

  % Parametros de Control de Iteracion
  tol_F = 1e-8;
  Maxiter = 10;

  i=1;     err_F = inf;
  
  % Iteracion de Metodo de Longitud de Arco.
  while and(i<Maxiter , err_F>tol_F)
    A = [F_x(x(k)+Dx_ki), -v ; 2*Dx_ki , 2*psi^2*v^2*DL_ki];
  
    b = [ -(f(x(k)+Dx_ki)-(L(k)+DL_ki)*v) ; -(Dx_ki^2+(psi*v*DL_ki)^2-Dl^2)];
  
    dIncr = A\b;    % solucion del sistema lineal
    
    Dx_ki = Dx_ki + dIncr(1); % Actualizo el Incremento en x
    DL_ki = DL_ki + dIncr(2); % Actualizo el Incremento en lambda
  
    err_F = norm(f(x(k)+Dx_ki)-(L(k)+DL_ki)*v); % error para control de convergencia
    
    i = i+1; % incrementa contador de iteraciones
  end

  % Tenemos incrementos con convergencia deseada para el paso(k)
  % Incremento x(k) para obtener x(k+1)
  x(k+1) = x(k) + Dx_ki;       L(k+1) = L(k) + DL_ki;  
end

ms=6; lw=2.2;

% Genera graficas de Solucion Numerica y Exacta
gridf = 0:.01:3;
Solexact=f(gridf);
plot(gridf,Solexact,'-r','markersize',ms,'linewidth',lw)
hold on
plot(x,L,'ok','markersize',ms,'linewidth',lw)
labx=xlabel('$x$'), laby=ylabel('$\lambda$')
%~ title('Solucion Con Metodo de Longitud de Arco')
hl =legend('Solución Exacta','Solución Numérica','location','North')
axis([0 3 -0.5 2.5],'square')%'equal')

set(gca, 'fontsize', 15 )

%~ copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(hl, "FontSize", 12);
set(labx, "FontSize", 14);
set(laby, "FontSize", 14);

print('Fig6Epscargabaja','-depslatex')
