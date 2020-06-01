##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
## _______     ___       __   __   __  |
##|   ____|   /   \     |  | |  | |  | | CATEDRA ESTRUCTURAS AERONAUTICAS III
##|  |__     /  ^  \    |  | |  | |  | |
##|   __|   /  /_\  \   |  | |  | |  | |    METODO DE RESIDUOS PONDERADOS
##|  |____ /  _____  \  |  | |  | |  | |
##|_______/__/     \__\ |__| |__| |__| | symTr: transpuesta variable symbolica
##                                     |
##---------CICLO LECTIVO 2019----------------------------------------------------

function d=symTr(z)

  d=sym(zeros(size(z,2),size(z,1)));
  
  for i=1:size(z,1)

    for j=1:size(z,2)

      d(j,i)=z(i,j);

    endfor
    
  endfor
  
endfunction
