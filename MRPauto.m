##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
## _______     ___       __   __   __  |
##|   ____|   /   \     |  | |  | |  | | CATEDRA ESTRUCTURAS AERONAUTICAS III
##|  |__     /  ^  \    |  | |  | |  | |
##|   __|   /  /_\  \   |  | |  | |  | |    METODO DE RESIDUOS PONDERADOS
##|  |____ /  _____  \  |  | |  | |  | |
##|_______/__/     \__\ |__| |__| |__| | MRPauto: script principal
##                                     |
##---------CICLO LECTIVO 2019----------------------------------------------------


## CONFIGURACION

pkg load symbolic; # Carga de paquete que me permite hacer operaciones algebraicas

warning ('off','OctSymPy:sym:rationalapprox');
## DECLARACIONES

syms x;

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

######### DOMINIO #############################

xInicial=-4;

xFinal=0;

######### OPERADOR DIFERENCIAL ################

function d=D(z)
  
  syms x;

  d=4*diff(z,x,2)-6*cos(x)*diff(z,x,2)+(x^2)*z;

endfunction

######## CONDICIONES HOMOGENEAS ###############

function d=forN(grados)
 
  d=sym(zeros(1,grados));
  
  for i=1:grados

    syms x;
    
    d(i)=(x+4)^(i+1); # Aca se carga el polinomio de interpolacion
       
  endfor
  
endfunction


######### CONDICIONES NO HOMOGENEAS ###########

psi=-5*x-17; # Cumple Dirichlet y Neumann

######### FUNCION F(X) ########################

fx=-4;

##//////// CONFIGURACION //////////////////////////////

gl=2 ;# Grados de libertad iniciales O subdominios

gli=gl;

maxIteraciones=10;

CRITERIO=15; # Criterio de convergencia en %

##//////// DECLARACIONES //////////////////////////////

CONVG=0;

N=sym(ones(1,1));

abscisas=linspace(xInicial,xFinal,10);

#############################################################
######### PUNTOS #######################################
#############################################################




figure(1);clf;

subplot(2,2,1)

hold on;

grid on

title ('COLOCACION POR PUNTOS')

while (CONVG<=0 && gl<=maxIteraciones)

######### GRADOS DE LIBERTAD ##################

  GL=linspace(xInicial,xFinal,gl);
  
  N=forN(gl);
    
######### VECTORES y MATRICES ############################
  
  F=subs(fx,x,GL);

  Fpsi=subs(D(psi),x,GL);

  DN=D(N);
  
  for i=1:gl
    
    try
      
      K=[K;subs(DN,x,GL(i))];

    catch

      K=[subs(DN,x,GL(i))];

    end_try_catch
  endfor
  
  
######### RESOLUCION ######################################

  C=inv(double(K))*(double(F)'-double(Fpsi)');

######### PLOTEO #########################################

    
  fp=function_handle(N*C+psi);
  
  AreaActual=trapz(abscisas,fp(abscisas));

  try
  
    AreaPrevia=trapz(abscisas,puntos);

  catch

    AreaPrevia=2*AreaActual;

  end_try_catch
  
  Error=abs((AreaPrevia-AreaActual)/AreaPrevia)*100;

  
  if  (Error<=CRITERIO || gl>=maxIteraciones)

    CONVG=1;

  endif
  
  puntos=fp(abscisas);

  plot(abscisas,puntos,["--" markStyle(gl) color(gl) ";GL= " num2str(gl) ";"]);
    
  gl++;
  
endwhile

hold off;

printf("\n\n\n COLOCACION POR PUNTOS %d GL -->  ERROR RELATIVO= %d%%\n",(gl-1),Error);



puntosSymbolico=N*C+psi

clear K Fpsi F;

N=sym(ones(1,1));

CONVG=0;

gp=gl;

gl=gli;


#############################################################
######### SUBDOMINIOS #######################################
#############################################################

subplot(2,2,2)

hold on;

grid on

title ('COLOCACION POR SUBDOMINIOS')

while (CONVG<=0 && gl<=maxIteraciones)
  
  GL=linspace(xInicial,xFinal,gl+1);

  N=forN(gl);
    
######### VECTORES y MATRICES ############################
  
  DPSI=D(psi);
  
  DN=D(N);
  
  for i=2:size(GL,2)
       
    try
      
      F=[F,integrador(fx,GL(i-1),GL(i))];

      Fpsi=[Fpsi,integrador(DPSI,GL(i-1),GL(i))];
      
      K=[K;integrador(DN,GL(i-1),GL(i))];

	
    catch
     
      F=[integrador(fx,GL(i-1),GL(i))];
      
      Fpsi=[integrador(DPSI,GL(i-1),GL(i))];

      K=[integrador(DN,GL(i-1),GL(i))];

      
    end_try_catch
  endfor
  
  
######### RESOLUCION ######################################

  C=inv(double(K))*(double(F)'-double(Fpsi)');

######### PLOTEO #########################################

  
  fp=function_handle(N*C+psi);
  
  AreaActual=trapz(abscisas,fp(abscisas));

  try
  
    AreaPrevia=trapz(abscisas,subdominios);

  catch

    AreaPrevia=2*AreaActual;

  end_try_catch
  
  Error=abs((AreaPrevia-AreaActual)/AreaPrevia)*100;
  
  if  (Error<=CRITERIO || gl>=maxIteraciones)

    CONVG=1;

  endif
  
  subdominios=fp(abscisas);
  
  plot(abscisas,subdominios,["--" markStyle(gl) color(gl) ";GL= " num2str(gl) ";"]);

  clear K Fpsi F;
  
  gl++;
  
endwhile

hold off;

printf("\n\n\n COLOCACION POR SUBDOMINIOS %d GL -->  ERROR RELATIVO= %d%%\n",(gl-1),Error);



subdominiosSymbolico=N*C+psi

clear K Fpsi F;

N=sym(ones(1,1));

CONVG=0;

gs=gl;

gl=gli;


#############################################################
######### GALERKIN #######################################
#############################################################

subplot(2,2,3)

hold on;

grid on

title ('GALERKIN')

while (CONVG<=0 && gl<=maxIteraciones)
  
  
  N=forN(gl);

  NT=symTr(N);

    
######### VECTORES y MATRICES ############################

  FX=NT*fx;
  
  DPSI=NT*D(psi);
  
  DN=NT*D(N);
  
  F=[integrador(FX,xInicial,xFinal)];

  Fpsi=[integrador(DPSI,xInicial,xFinal)];

  K=[integrador(DN,xInicial,xFinal)];

  
  
######### RESOLUCION ######################################

  C=inv(double(K))*(double(F)-double(Fpsi));

######### PLOTEO #########################################

  fp=function_handle(N*C+psi);
  
  AreaActual=trapz(abscisas,fp(abscisas));

  try
  
    AreaPrevia=trapz(abscisas,galerkin);

  catch

    AreaPrevia=2*AreaActual;

  end_try_catch
  
  Error=abs((AreaPrevia-AreaActual)/AreaPrevia)*100;
  
  if  (Error<=CRITERIO || gl>=maxIteraciones)

    CONVG=1;

  endif

  RESULTADO=N*C+psi;
  
  galerkin=fp(abscisas);
  
  plot(abscisas,galerkin,["--" markStyle(gl) color(gl) ";GL= " num2str(gl) ";"]);

  clear K Fpsi F NT;

  N=sym(ones(1,1));

  gl++;
  
endwhile

hold off;

printf("\n\n\n GALERKIN %d GL -->  ERROR RELATIVO= %d%%\n",(gl-1),Error);



galerkinSymbolico=RESULTADO

CONVG=0;

gg=gl;

gl=gli;

#############################################################
######### CUADRADOS MINIMOS #################################
#############################################################

subplot(2,2,4)

hold on;

grid on

title ('CUADRADOS MINIMOS')

while (CONVG<=0 && gl<=maxIteraciones)
  

  N=forN(gl);

  NT=symTr(N);
  
    
######### VECTORES y MATRICES ############################

  DNT=D(NT);
  
  FX=DNT*fx;
  
  DPSI=DNT*D(psi);
  
  DN=DNT*D(N);
  
  F=[integrador(FX,xInicial,xFinal)];

  Fpsi=[integrador(DPSI,xInicial,xFinal)];

  K=[integrador(DN,xInicial,xFinal)];

  
  
######### RESOLUCION ######################################

  C=inv(double(K))*(double(F)-double(Fpsi));

######### PLOTEO #########################################

  fp=function_handle(N*C+psi);
  
  AreaActual=trapz(abscisas,fp(abscisas));

  try
  
    AreaPrevia=trapz(abscisas,cuadrados);

  catch

    AreaPrevia=2*AreaActual;

  end_try_catch
  
  Error=abs((AreaPrevia-AreaActual)/AreaPrevia)*100;
  
  if  (Error<=CRITERIO || gl>=maxIteraciones)

    CONVG=1;

  endif

  RESULTADO=N*C+psi;
  
  cuadrados=fp(abscisas);
  
  plot(abscisas,cuadrados,["--" markStyle(gl) color(gl) ";GL= " num2str(gl) ";"]);

  clear K Fpsi F NT;
  
  N=sym(ones(1,1));

  gl++;
  
endwhile

hold off;

printf("\n\n\n CUADRADOS MINIMOS %d GL -->  ERROR RELATIVO= %d%%\n",(gl-1),Error);



cuadradosSymbolico=RESULTADO

CONVG=0;

############# SOLUCIONES ###################

figure (2);clf;

hold on;

grid on

title ('SOLUCIONES')

plot(abscisas,puntos,["--" markStyle(1) color(1) ";PUNTOS - GL= " num2str(gp-1) ";"]);

plot(abscisas,subdominios,["--" markStyle(2) color(2) ";SUBDOMINIOS - GL= " num2str(gs-1) ";"]);

plot(abscisas,galerkin,["--" markStyle(3) color(3) ";GALERKIN - GL= " num2str(gg-1) ";"]);

plot(abscisas,cuadrados,["--" markStyle(4) color(4) ";CUADRADOS - GL= " num2str(gl-1) ";"]);

hold off
