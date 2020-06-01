##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
## _______     ___       __   __   __  |
##|   ____|   /   \     |  | |  | |  | | CATEDRA ESTRUCTURAS AERONAUTICAS III
##|  |__     /  ^  \    |  | |  | |  | |
##|   __|   /  /_\  \   |  | |  | |  | |    METODO DE RESIDUOS PONDERADOS
##|  |____ /  _____  \  |  | |  | |  | |
##|_______/__/     \__\ |__| |__| |__| | integrador: integrac. metodo trapezoides
##                                     |
##---------CICLO LECTIVO 2019----------------------------------------------------


function d=integrador(z,xi,xf)
  
  PRESICION=10;
  
  dom=linspace(xi,xf,PRESICION);
  
  d=zeros(size(z));
  
  for i=1:size(z,1)

    for j=1:size(z,2)	     

      test=sym(z(i,j));
      
      if (isempty(symvar(test))==1)

	y=ones(size(dom))*double(z(i,j)); # Si el integrando es constante

      else
      
	fp=function_handle(z(i,j));
	
	y=fp(dom);
	
      endif  
	
	
      try
	
	d(i,j)=trapz(dom,y);

      catch  # Si es nulo

	d(i,j)=0;

      end_try_catch

      
    endfor
    
  endfor
  

endfunction

