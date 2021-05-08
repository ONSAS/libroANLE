clear all, close all

% parametros de geometria materiales y secciones de Seccion de apuntes.
Es =[ 5e5 ]; As = [ 100 ]; auxy = 2500 ; auxx = 2500; Pext = -2*3.0e6 ;
% parametros de geometria materiales y secciones de Ejemplo de articulo de Li.
%~ Es=[ 210e9 ]; As=[ 2.5e-4 ]; auxx =cos(65*pi/180)*2; auxy =sin(65*pi/180)*2; Pext= -1.0e6;

% conectividad
Nodes = [ 0 0 ; ...
          auxx auxy ; ...
          2*auxx 0 ] ;
l0 = sqrt(auxx^2+auxy^2) ;

Conec = [ 1 2 1 1 ;
          2 3 1 1 ] ;

% calcula numero de nodos y elementos
nnodes = size(Nodes,1); nelems = size(Conec,1);

% grados de libertad fijos
fixeddofs = [ 1 2 5 6 ];

% carga externa
Fext    = zeros(2*nnodes,1) ; Fext(4) = Pext ;

% tolerancias criterio de parada
tolk  = 50;  tolu = 1e-4;

% --- calculos previos -----------------------
Ge = [ 1 0 -1 0; 0 1 0 -1; -1 0 1 0; 0 -1 0 1 ] ;

% calcula los grados de libertad libres
freedofs = (1:(2*nnodes))'; freedofs(fixeddofs) = [];

% calcula largos, angulos y cosenos de angulos que forman barras con x
largosini = sqrt( ( Nodes( Conec(:,2),1) - Nodes( Conec(:,1),1) ).^2 ...
                + ( Nodes( Conec(:,2),2) - Nodes( Conec(:,1),2) ).^2 ) ;
thetasini = atan2( ( Nodes( Conec(:,2),2) - Nodes( Conec(:,1),2) ) , ...
                   ( Nodes( Conec(:,2),1) - Nodes( Conec(:,1),1) ) ) ;
cosini   = cos( thetasini ) ; sinini  = sin( thetasini ) ;
xelems   = reshape( Nodes( Conec(:,1:2) ,1 )', nelems,2 ) ;
yelems   = reshape( Nodes( Conec(:,1:2) ,2 )', nelems,2 ) ;

% --- iteracion de newton raphson -------
% inicializa Uk con un vector de zeros
Uk      = zeros(2*nnodes,1);   Stressk = zeros(  nelems,1);  
Fext(fixeddofs ) = [] ; histuks = [];

fin = 0 ; k = 0;
fprintf('iter &    $u^k(4)$  &$\\varepsilon^k$&  $\\sigma^k$ & fin? \\\\ \n \\hline')
while fin == 0,
  k += 1; % suma 1 al contador de iteraciones
  KTk   = sparse( 2*nnodes , 2*nnodes ) ;  Fintk = zeros ( 2*nnodes , 1        ) ;
  
  for elem = 1:nelems
    nodeselem = Conec(elem,1:2)' ;   dofselem  = nodes2dofs( nodeselem , 2 ) ;    
    E  = Es(Conec(elem,3)) ;  A  = As(Conec(elem,4)) ;  le = largosini(elem)   ;
  
    Xe = [ xelems( elem, 1) yelems( elem, 1) xelems( elem, 2) yelems( elem, 2) ]' ;
    Ue = Uk(dofselem) ;
    
    % calcula vectores b1 y b2
    B1e = 1.0 / ( le^2 ) * Xe' * Ge ;  B2e = 1.0 / ( le^2 ) * Ue' * Ge ;
    
    % calcula matrices de rigidez y fuerzas internas
    KT1e  = E * A * le * ( B1e' * B1e ) ;
    KT2e  = E * A * le * ( B2e' * B1e + B1e' * B2e + B2e' * B2e ) ;
    Ksige = A * Stressk(elem) / le * Ge  ;
    Finte = (B1e+B2e)' * A * le * Stressk(elem) ;
    
    % ensamblado de matrices y vector de fuerzas internas
    KTk   ( dofselem, dofselem) = KTk (dofselem,dofselem) + KT1e + KT2e + Ksige       ;
    Fintk ( dofselem          ) = Fintk(dofselem) + Finte ;
  end

  KTk(fixeddofs, :) = [] ;  KTk(:, fixeddofs) = [] ;  Fintk(fixeddofs ) = [] ;

  DeltaUk      = KTk \ ( Fext - Fintk )  ;
  normUk       = norm( Uk(freedofs) ) ;
  Uk(freedofs) = Uk(freedofs) + DeltaUk ;
  
  for elem = 1:nelems
    nodeselem = Conec(elem,1:2)' ;  dofselem  = nodes2dofs( nodeselem , 2 ) ;
    E  = Es(Conec(elem,3)) ;  A  = As(Conec(elem,4)) ;  le = largosini(elem)   ;

    Xe = [ xelems( elem, 1) yelems( elem, 1) xelems( elem, 2) yelems( elem, 2) ]' ;
    Ue = Uk(dofselem) ;
    
    B1e = 1.0 / ( le^2 ) * Xe' * Ge ;  B2e = 1.0 / ( le^2 ) * Ue' * Ge ;

    epsGelem = ( B1e + 0.5*B2e ) * Ue ;    Stressk(elem) = E * epsGelem ;
  end

  if ( k > tolk ) || ( norm( DeltaUk )  < tolu * normUk ), fin = 1; end

  histuks = [ histuks; Uk(4) ] ;
  fprintf(' %3i & %12.3e & %12.3e & %12.3e & %1i \\\\ \n', k,  Uk(4), Stressk(1)/E, Stressk(1), fin )
end

NodesDef = Nodes + reshape( Uk , 2,nnodes)' ;
xelemsdef = reshape( NodesDef( Conec(:,1:2) ,1 )', nelems,2 ) ;
yelemsdef = reshape( NodesDef( Conec(:,1:2) ,2 )', nelems,2 ) ;

w = Uk(4) ; PGTeo = 2*Es*As*(auxy+w)*(2*auxy*w+w^2) / ( 2*l0^3) 
            PGNum = 2*Stressk(1)*A/l0*(auxy+w)

Error_relativo_PGNum = abs( Pext - PGNum ) / abs(Pext)
Error_relativo_PGTeo = abs( Pext - PGTeo ) / abs(Pext)

plotfontsize = 22 ;

% --- ploteos ---
figure, hold on
plot(xelems,    yelems,    'b--o','linewidth',4.5,'markersize',13)
plot(xelemsdef, yelemsdef, 'r-x','linewidth',4.5,'markersize',13)
labx=xlabel('x'), laby=ylabel('z'), axis equal

set(gca, 'fontsize', plotfontsize )
set(labx, "FontSize", plotfontsize); set(laby, "FontSize", plotfontsize);
print('figmises','-depslatex')

figure
plot(-histuks, 'b-x','linewidth',4.5,'markersize',13)
labx=xlabel('Iteraciones'), laby=ylabel('|$u^k(4)$|')
set(gca, 'fontsize', plotfontsize )
set(labx, "FontSize", plotfontsize); set(laby, "FontSize", plotfontsize);
print('histukmises','-depslatex')
