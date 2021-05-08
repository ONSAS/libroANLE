%% Medidas de Deformacion No-Lineales %% Grafica de Resultados de Seccion 2.3
clc, clear all, close all

% Parametros de la Estructura
x  = 2500; z = 2500;
A0 = 100 ; E = 5e5 ;

% Calculo longitud incial barra
l0 = sqrt(z^2+x^2);

% Se define q para las distintas medidas de deformacion
PE = @(w) E*A0*(z+w).*(sqrt((w+z).^2+x^2)-l0)./(l0*sqrt((w+z).^2+x^2));
PG = @(w) E*A0*(z+w).*(2*z*w+w.^2)/(2*l0^3);
PL = @(w) E*A0*l0*(z+w).*log(sqrt((z+w).^2+x^2)/l0)./((z+w).^2+x^2);

% Se define rango de valores de w para graficar
w = linspace(-6000,0,100);  wref1=w(40); wref2=w(80); wref3=w(1);

% Se grafica
ms=12; lw=2.5;
figure, hold on
plot(-wref1,-PE(wref1)/1e6,'s-r','markersize',ms,'linewidth',lw)
plot(-wref1,-PG(wref1)/1e6,'*-b','markersize',ms,'linewidth',lw)
plot(-wref1,-PL(wref1)/1e6,'o-k','markersize',ms,'linewidth',lw)

plot(-wref2,-PE(wref2)/1e6,'s-r','markersize',ms,'linewidth',lw)
plot(-wref2,-PG(wref2)/1e6,'*-b','markersize',ms,'linewidth',lw)
plot(-wref2,-PL(wref2)/1e6,'o-k','markersize',ms,'linewidth',lw)

plot(-wref3,-PE(wref3)/1e6,'s-r','markersize',ms,'linewidth',lw)
plot(-wref3,-PG(wref3)/1e6,'*-b','markersize',ms,'linewidth',lw)
plot(-wref3,-PL(wref3)/1e6,'o-k','markersize',ms,'linewidth',lw)

plot(-w,-PE(w)/1e6,'-r','markersize',ms,'linewidth',lw)
plot(-w,-PG(w)/1e6,'-b','markersize',ms,'linewidth',lw)
plot(-w,-PL(w)/1e6,'-k','markersize',ms,'linewidth',lw)
legend('$P_E$: Def. Unit. Ing. Rotada','$P_G$: Def. Unit. de Green','$P_L$: Def. Unit. Log. Rotada','location','NorthWest')
labx=xlabel('Desplazamiento: (-$w$)'); laby=ylabel('Carga: (-$P$ $\times \, 10^{+6}$)')
%~ set(gca, 'linewidth', 2, 'fontsize', 15 )
set(gca, 'fontsize', 15 )

copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 13);
set(labx, "FontSize", 14);
set(laby, "FontSize", 14);

print( 'nonlindef','-depslatex') ;

