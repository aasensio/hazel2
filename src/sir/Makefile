# ****************************************************************************
# 
# Fichero   : Makefile
# Proyecto  : Metodos de inversion
# Dpto.     : Departamento de Astrofisica - ULL
# Autor     : Basilio Ruiz Cobo - brc@ll.iac.es
# Proposito : Mantenimiento de los ejecutables de inversion
#
# Uso       : > make [fc=compilador]               ... genera sir.x
#             > make [fc=compilador] lande.x       ... genera lande.x
#             > make [fc=compilador] modelador3.x  ... genera modelador3.x
#
#             donde fc puede ser gfortran, ifort o nada (defecto = ifort).
#             Por ejemplo:
#
#             > make sir.x
#             > make fc=gfortran lande.x
#             > make fc=ifort modelador3.x
#
# 27/01/04  : primera version.
# 12/11/15  : el compilador a usar puede entrar por linea de comandos.
#
# ****************************************************************************


# Recordatorio de la ultima compilacion.

.KEEP_STATE:

# Compilador, flags del compilador y flags de construccion de la libreria.
#
# Extension de lineas de codigo: -ffixed-line-length-none o -extend_source

FC = gfortran
FFLAGS = -O3 -ffixed-line-length-none -fno-automatic -fPIC
ARFLAGS = ruv
 
# Fichero de configuracion (definicion de dimensiones).

FILECONF = PARAMETER
FILEMAKE = .make.lib

# Ficheros fuentes de los ejecutables.

MAIN1SRC = StokesFR.f
#MAIN2SRC = modelador2.f
#MAIN3SRC = lande.f
#MAIN4SRC = modelador3.f

# Libreria de inversion y todos sus modulos.

LIB = libreria.a

LIBSRC =atmdatb.f \
	cabecera.f \
	chachi.f \
	conhsra.f \
	deconv.f \
	deconv2.f \
	delta1.f \
	deriva.f \
	deriva4cuad.f \
	derivacuad.f \
	dsaha.f \
	dtplanck.f \
	equi.f \
	escribe_stokes.f \
	fint.f \
	gasb.f \
	gase.f \
	hermite.f \
	hermite_c.f \
	intmodel.f \
	kappach.f \
	leeabun.f \
	leelineasii.f \
	lee_malla.f \
	lee_model.f \
	leyendo.f \
	locate.f \
	lubksb.f \
	ludcmp.f \
	matinx.f \
	matinx2.f \
	mensaje.f \
	meves.f \
	molecb.f \
	mreadi2.f \
	mreadi3.f \
	mreadi33.f \
	mvoigt.f \
	mvoigtc.f \
	planck.f \
	plank2.f \
	polint.f \
	refrax.f \
	rnorma.f \
	saha.f \
	sublande.f \
	sublande2.f \
	StokesFRsub.f \
	taulinea.f \
	voigt.f \
	zeeman.f \
	zeeman_jk.f \
	interface.f

# Referencias para el reemplazamiento de sufijos.

MAIN1OBJ = ${MAIN1SRC:.f=.o}
#MAIN2OBJ = ${MAIN2SRC:.f=.o}
#MAIN3OBJ = ${MAIN3SRC:.f=.o}
#MAIN4OBJ = ${MAIN4SRC:.f=.o}

LIBOBJ = ${LIBSRC:.f=.o,.f90=.o}


#
# A)
# Ejecutable sir.x ('uso> make' o 'uso> make sir.x').
#

StokesFR.x: ${MAIN1OBJ} ${LIB} ${FILEMAKE}
	${LINK.f} -o $@ ${MAIN1OBJ} ${LIB}

clean: 
	rm -rf *.o *.a


#
# B)
# Ejecutable modelador2.x ('uso> make modelador2.x').
#
#
#modelador2.x: ${MAIN2OBJ} ${LIB}
#	${LINK.f} -o $@ ${MAIN2OBJ} ${LIB}

#
# C)
# Ejecutable lande.x ('uso> make lande.x').
#
#
#lande.x: ${MAIN3OBJ} ${LIB}
#	${LINK.f} -o $@ ${MAIN3OBJ} ${LIB}

# B2)
# Ejecutable modelador3.x ('uso> make modelador3.x').
#
#
#modelador3.x: ${MAIN4OBJ} ${LIB}
#	${LINK.f} -o $@ ${MAIN4OBJ} ${LIB} 


#
# D)
# Actualizacion de la libreria.
#

${LIB}: ${LIBOBJ}
	${AR} ${ARFLAGS} $@ $?


#
# E)
# Reconstruccion de la libreria.
#

${FILEMAKE}: ${FILECONF}
	${COMPILE.f} ${LIBSRC}
	${AR} ${ARFLAGS} ${LIB} ${LIBOBJ}
	@touch $@


# ****************************************************************************
