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
        !fichabun='THEVENIN'        

        !call leeabun(0,eps)
        eps(1) = 12.00
        eps(2) = 11.00
        eps(3) = 1.00
        eps(4) = 1.15
        eps(5) = 2.60
        eps(6) = 8.69
        eps(7) = 7.99
        eps(8) = 8.91
        eps(9) = 4.56 
        eps(10) = 8.00 
        eps(11) = 6.28 
        eps(12) = 7.53 
        eps(13) = 6.43 
        eps(14) = 7.50 
        eps(15) = 5.45 
        eps(16) = 7.21 
        eps(17) = 5.50 
        eps(18) = 6.58 
        eps(19) = 5.05 
        eps(20) = 6.36 
        eps(21) = 2.99 
        eps(22) = 4.88 
        eps(23) = 3.91 
        eps(24) = 5.61 
        eps(25) = 5.47 
        eps(26) = 7.46 
        eps(27) = 4.85 
        eps(28) = 6.18 
        eps(29) = 4.24 
        eps(30) = 4.60 
        eps(31) = 2.88 
        eps(32) = 3.57 
        eps(33) = 2.39 
        eps(34) = 3.35 
        eps(35) = 2.63 
        eps(36) = 3.21 
        eps(37) = 2.60 
        eps(38) = 2.93 
        eps(39) = 2.18 
        eps(40) = 2.46 
        eps(41) = 1.46 
        eps(42) = 2.10 
        eps(43) = 0.00 
        eps(44) = 1.78 
        eps(45) = 1.10 
        eps(46) = 1.69 
        eps(47) = 0.94 
        eps(48) = 1.86 
        eps(49) = 1.66 
        eps(50) = 2.00 
        eps(51) = 1.00 
        eps(52) = 2.25 
        eps(53) = 1.51 
        eps(54) = 2.19 
        eps(55) = 1.12 
        eps(56) = 2.18 
        eps(57) = 1.07 
        eps(58) = 1.58 
        eps(59) = 0.76 
        eps(60) = 1.40 
        eps(61) = 0.00 
        eps(62) = 0.88 
        eps(63) = 0.48 
        eps(64) = 1.13 
        eps(65) = 0.20 
        eps(66) = 1.07 
        eps(67) = 0.26 
        eps(68) = 0.93 
        eps(69) = 0.00 
        eps(70) = 1.08 
        eps(71) = 0.76 
        eps(72) = 0.88 
        eps(73) = -.09 
        eps(74) = 0.98 
        eps(75) = 0.26 
        eps(76) = 1.45 
        eps(77) = 1.36 
        eps(78) = 1.80 
        eps(79) = 1.13 
        eps(80) = 1.27 
        eps(81) = 0.90 
        eps(82) = 1.90 
        eps(83) = 0.71 
        eps(84) = -8.0 
        eps(85) = -8.0 
        eps(86) = -8.0 
        eps(87) = -8.0 
        eps(88) = -8.0 
        eps(89) = -8.0 
        eps(90) = 0.02 
        eps(91) = -8.0 
        eps(92) = -.47
        
c leemos la malla  
        mallaobs='malla.grid'
        call lee_malla(mallaobs)     ! carga el common Malla/ntl,nlin,npas,nble,dlamda
                          
c nombre del fichero de salida con los perfiles de Stokes        
        Stokesfilename='perfil.per'        
        
        return      
        end