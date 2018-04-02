        subroutine leyendo !carga en commons toda la informacion necesaria
               
c para definir los datos de entrada
        character*100 mallaobs
        character*100 nomlineas,fichabun,modelname
        character*100 Stokesfilename
        logical ex        
        real*4 eps(92)
                       
c lugares comunes de memoria para 
        common/ficlineas/nomlineas           !para leelineasii (en StokesFRsub)
        common/fichabun/fichabun
        common/OutputStokes/Stokesfilename   !para escribe_Stokes
        common/abundances/eps
        
c nombre del fichero con los parametros atomicos
        nomlineas='LINEAS'
        
c nombre del fichero con las abundancias
        fichabun='THEVENIN'        

        call leeabun(0,eps)
        
c leemos la malla  
        mallaobs='malla.grid'
        call lee_malla(mallaobs)     ! carga el common Malla/ntl,nlin,npas,nble,dlamda
                          
c nombre del fichero de salida con los perfiles de Stokes        
        Stokesfilename='perfil.per'        
        
        return      
        end