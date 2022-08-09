% Ejemplo de Metodo de Diferencia Centrada: Edificio Bajo Carga de Explosion %
clc, clear

% Definicion de la estructura
A     = 30; %m				Lado del Edificio
esp   = 0.25; %m			Espesor Losa Piso
hv    = 0.7; %m			Altura Viga
bv    = 0.3; %m			Ancho Viga
muro  = 300; %kg/m2		Peso muro por uni. area
ncols = 16; %			Cantidad de Columnas (4x4)

a = 0.7; %m				Lado de Columna
E = 30e9; %N/m2			Mod. Young Hormigon
I = a^4/12; %m4
H = 3.5; %m				Altura Entre Pisos

% Definicion de Matriz de Rigidez
k = ncols * 12*E*I/H^3; %N/m
       %u1 u2 u3
K = k*[2 -1 0 ; -1 2 -1 ; 0 -1 1]; %N/m

% Definicion de Matriz de Masa
Mpiso  = 2500 * ( esp*A^2+hv*bv*A*8 ) ; %kg
Mpilar = 2500*ncols*a^2*H; %kg
Mmuro  = muro*H*A*4*0.6; %kg
Mint   = Mpiso + Mpilar + Mmuro;
Msup   = Mpiso + Mpilar/2 + Mmuro/2;

M = diag([Mint, Mint, Msup]); %kg

% Definicion de Matriz Amortiguamiento
    % Usamos Rayleigh Damping: C = eta*M + delta*K
    % Amortiguamiento: 3% critico para 25rad/s y 106rad/s
C = 1.21*M + 4.6e-4*K;

% Calculo de Modos Normales
[PHI, w2] = eig(K,M);

w = sqrt(diag(w2)); %rad/s
f = w/2/pi; %Hz
T = 1./f %sec

% Calculo de Paso Critico para Diferencia Centrada
dtcr = min(T)/pi

% Definicion de Historia de Presiones y Fuerzas Laterales
Aint = A*H; %m2
Asup = A*H/2; %m2

% funcion analitica para evaluar presion triangular en el tiempo
press = @(t,ta,te,pr1,pr2) (pr1-(pr1-pr2)/te*(t-ta)).*and(t>ta,t<ta+te);

% historia de presion en el tiempo como suma de tramos lineales (triangulos)
pt = @(t) press(t,15.6e-3,6.2e-3,0.55e6,0)+ ...
          press(t,66.0e-3,11.9e-3,0,-0.035e6)+...
          press(t,77.9e-3,19.1e-3,-0.035e6,0); %N/m2

ft = @(t) pt(t)*[Aint;Aint;Asup]; %N

% Definicion de Condiciones Iniciales y Tiempo Final
t0 = 0; %sec
u0 = [0;0;0]; %m
v0 = [0;0;0]; %m
ac0 = M\(ft(t0)-C*v0-K*u0); % de ec de movimiento Mu.. + Ku = ft

tf = 2*max(T); % 2 periodos del modo fundamental

% Inicializacion Difrerencia Centrada
dt = dtcr/20;
a0 = 1/dt^2; a1=1/2/dt; a2=2*a0; a3=1/a2;

u(:,1) = u0 - dt*v0 + a3*ac0; % u(-dt)
u(:,2) = u0; % u(0)

Meff = a0*M+a1*C;   Keff = (K-a2*M);   M2 = a0*M-a1*C;

% Comienza Marcha en el Tiempo usando Diferencia Centrada
t(1) = t0-dt;   t(2) = t0;   k=2;
while t(end) < tf
    feff = ft(t(k)) - Keff*u(:,k) - M2*u(:,k-1);
    u(:,k+1) = Meff\feff;
    acc(:,k) = a0*(u(:,k+1)-2*u(:,k)+u(:,k-1));
    vel(:,k) = a1*(u(:,k+1)-u(:,k-1));
    k=k+1;
    t(k)=t(k-1)+dt;
end

subplot(1,2,1),    plot(t,pt(t)/1e3,'-b','linewidth',2)
labx=xlabel('t [s]'), laby=ylabel('presion [kPa]')
axis([t0,0.3*max(t),1.2*min(min(pt(t)/1e3)),1.2*max(max(pt(t)/1e3))])
%~ title('Presion por Explosion')
fontsize = 13;

set(labx, "FontSize", fontsize); set(laby, "FontSize",fontsize);
set(gca, 'fontsize', fontsize )

fontsize = 15;

subplot(1,2,2), hold on
plot(t,u(1,:),'-b','linewidth',2)
plot(t,u(2,:),'r--','linewidth',2)
plot(t,u(3,:),'-k','linewidth',2)
hl = legend('$u_1$','$u_2$','$u_3$'), labx=xlabel('t [s]'), laby=ylabel('desplazamiento [m]')
axis([t0,1.1*max(t),1.5*min(min(u)),1.5*max(max(u))])
%~ title('Respuesta de la Estructura')
set(hl,'fontsize',13);
set(labx, "FontSize", fontsize); set(laby, "FontSize", fontsize);
set(gca, 'fontsize', fontsize )

print('Explicit','-depslatex');
%print('Explicit','-dpdflatex');
