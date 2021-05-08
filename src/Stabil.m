% Codigo de ejemplo de implementacion de metodo de analisis elastoplastico
% Funcion de endurecimiento lineal. Historia de deformaciones dada.

z = 0:.01:3;

lam1 = (2-z.^2)/2 +sqrt( (2-z.^2).^2 / 4 -1 ) ;
lam2 = (2-z.^2)/2 -sqrt( (2-z.^2).^2 / 4 -1 ) ;

maxs = max( abs( [  lam1; lam2] ) ) ;

figure
lw = 3.5; ms = 5;
plot(z, maxs,'b','linewidth',lw)
labx=xlabel('$\omega_i \Delta t$'), laby=ylabel('$\rho$')
set(labx, "FontSize", 20); set(laby, "FontSize", 20);
set(gca, 'fontsize', 20 )
print('Stabil','-dpdflatex')
