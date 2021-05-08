%% Medidas de Deformacion No-Lineales %% Grafica de Resultados de Seccion 2.3
clc; clear

% Parametros de la Estructura
x=2500; z=2500;
A0 = 100; E = 5e5;

% Calculo longitud incial barra
l0 = sqrt(z^2+x^2);

% Se define q para las distintas medidas de deformacion
PE = @(w) E*A0*(z+w).*(sqrt((w+z).^2+x^2)-l0)./(l0*sqrt((w+z).^2+x^2));
PG = @(w) E*A0*(z+w).*(2*z*w+w.^2)/(2*l0^3);
PL = @(w) E*A0*l0*(z+w).*log(sqrt((z+w).^2+x^2)/l0)./((z+w).^2+x^2);

% Se define rango de valores de w para graficar
w = linspace(-6000,0,100);

% Se grafica
plot(-w,-PE(w)/1e6,'--r',-w,-PG(w)/1e6,'-b',-w,-PL(w)/1e6,':k')
legend('P_E: Def. Unit. Ing. Rotada','P_G: Def. Unit. de Green','P_L: Def. Unit. Log. Rotada')
xlabel('Desplazamiento: (-w)'); ylabel('Carga: (-P x 10^+^6)')
