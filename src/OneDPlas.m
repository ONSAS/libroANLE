% Codigo de ejemplo de implementacion de metodo de analisis elastoplastico
% Funcion de endurecimiento lineal. Historia de deformaciones dada.
clear all, close all

% parametros constitutivos
E      = 210e9 ; sigmaY = 250e6 ;
%~ K      =  0 ;
K      =  21e9 ;

% historia de deformaciones impuesta
epsmax       = sigmaY/E * 1.5 ; 
deps         = epsmax/20 ; % incremento 
loadeps      = (0:deps:epsmax)' ; % carga
perunloadeps = (epsmax:-deps:-epsmax)' ; % descarga
perloadeps   = (-epsmax:deps:epsmax)'  ; % carga
epshist      = [ loadeps; perunloadeps; perloadeps; perunloadeps ] ;
ntimes       = length(epshist);

% vectores de historia tension, deformacion plastica, elastica y acumulada
epsplhist   = zeros(ntimes,1) ;   epsplachist = zeros(ntimes,1) ; 
epselhist   = zeros(ntimes,1) ;   sigmahist   = zeros(ntimes,1) ;

% se considera que en el tiempo inicial (1) todas las magnitudes son nulas
for i=2:ntimes
  epsetrial  = epshist(i) - epsplhist(i-1)  ;
  sigmatrial = E * epsetrial  ;
  epsplactrial = epsplachist(i-1)  ;
  
  phitrial   = abs(sigmatrial) - ( sigmaY  + K*epsplactrial ) ;
  
  if phitrial <= 0 % se continua en rango elastico
    sigmahist(i)   = sigmatrial ;
    epsplhist(i)   = epsplhist(i-1); 
    epsplachist(i) = epsplachist(i-1);
    Cep            = E ;
  
  else % se debe calcular la deformacion plastica
    Deltagamma     = phitrial / (E + K) ;
    sigmahist(i)   = (1-Deltagamma*E/abs(sigmatrial))* sigmatrial  ;
    epsplhist(i)   = epsplhist(i-1) + Deltagamma * sign(sigmatrial) ; 
    epsplachist(i) = epsplachist(i-1) + Deltagamma ;
    Cep            = E*K/(E+K) ;
  end
end

figure, grid on, lw = 3; ms = 5;
plot(epshist, sigmahist,'b-x','linewidth',lw,'markersize',ms)
labx=xlabel('Deformacion'), laby=ylabel('Tension')
set(labx, "FontSize", 20); set(laby, "FontSize", 20);
set(gca, 'linewidth', 2, 'fontsize', 20 ), print('plasticcycle2','-depslatex')
