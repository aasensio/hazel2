project(
  'hazel2',
  'c', 'cython', 'fortran',
  # The single source of the package version. Bump it here, then tag vX.Y.Z.
  version: '1.1.0',
  license: 'MIT',
  meson_version: '>=1.4.0',
  default_options: ['buildtype=release'],
)

py = import('python').find_installation(pure: false)
numpy_dep = dependency('numpy')

# Binary wheels have to run on any CPU, so the default is a portable build.
# For a local build tuned to this machine:
#   pip install . -Csetup-args=-Dnative=true
native_args = get_option('native') ? ['-march=native'] : []

free_form_args = ['-ffree-form', '-ffree-line-length-none', '-w'] + native_args
fixed_form_args = ['-fno-automatic', '-ffixed-line-length-none', '-w'] + native_args

# Hazel
py.extension_module(
  'hazel_code',
  'src/hazel/vars.f90',
  'src/hazel/singleton.f90',
  'src/hazel/maths.f90',
  'src/hazel/allen.f90',
  'src/hazel/svd.f90',
  'src/hazel/io_py.f90',
  'src/hazel/SEE.f90',
  'src/hazel/rt_coef.f90',
  'src/hazel/synth.f90',
  'src/hazel/transfer.f90',
  'src/hazel/hazel_py.f90',
  'src/hazel/hazel_code.pyx',
  fortran_args: free_form_args,
  dependencies: numpy_dep,
  subdir: 'hazel/codes',
  install: true,
)

# SIR. The fixed-form sources live in their own library because Meson applies
# fortran_args per target, and -fno-automatic must not reach interface.f90.
sir_f77 = static_library(
  'sir_f77',
  'src/sir/StokesFR.f',
  'src/sir/StokesFRsub.f',
  'src/sir/atmdatb.f',
  'src/sir/cabecera.f',
  'src/sir/chachi.f',
  'src/sir/conhsra.f',
  'src/sir/deconv.f',
  'src/sir/deconv2.f',
  'src/sir/delta1.f',
  'src/sir/deriva.f',
  'src/sir/deriva4cuad.f',
  'src/sir/derivacuad.f',
  'src/sir/dsaha.f',
  'src/sir/dtplanck.f',
  'src/sir/equisubmu.f',
  'src/sir/escribe_stokes.f',
  'src/sir/fint.f',
  'src/sir/gasb.f',
  'src/sir/gase.f',
  'src/sir/hermite.f',
  'src/sir/hermite_c.f',
  'src/sir/inicia_pefrompgt.f',
  'src/sir/intmodel.f',
  'src/sir/kappach.f',
  'src/sir/lee_all_lines.f',
  'src/sir/lee_malla.f',
  'src/sir/lee_model.f',
  'src/sir/leeabun.f',
  'src/sir/leelineasii.f',
  'src/sir/leyendo.f',
  'src/sir/locate.f',
  'src/sir/lubksb.f',
  'src/sir/ludcmp.f',
  'src/sir/matinx.f',
  'src/sir/matinx2.f',
  'src/sir/mensaje.f',
  'src/sir/meves.f',
  'src/sir/molecb.f',
  'src/sir/mreadi2.f',
  'src/sir/mreadi3.f',
  'src/sir/mreadi33.f',
  'src/sir/mvoigt.f',
  'src/sir/mvoigtc.f',
  'src/sir/planck.f',
  'src/sir/plank2.f',
  'src/sir/polint.f',
  'src/sir/refrax.f',
  'src/sir/rnorma.f',
  'src/sir/saha.f',
  'src/sir/sublande.f',
  'src/sir/sublande2.f',
  'src/sir/taulinea.f',
  'src/sir/voigt.f',
  'src/sir/zeeman.f',
  'src/sir/zeeman_jk.f',
  fortran_args: fixed_form_args,
  include_directories: include_directories('src/sir'),
  pic: true,
)

py.extension_module(
  'sir_code',
  'src/sir/interface.f90',
  'src/sir/sir_code.pyx',
  fortran_args: free_form_args,
  include_directories: include_directories('src/sir'),
  link_whole: sir_f77,
  dependencies: numpy_dep,
  subdir: 'hazel/codes',
  install: true,
)

# Python package
py.install_sources(
  'hazel/__init__.py',
  'hazel/analysis.py',
  'hazel/atmosphere.py',
  'hazel/chromosphere.py',
  'hazel/configuration.py',
  'hazel/exceptions.py',
  'hazel/forward_nn.py',
  'hazel/forward_nn_transformer.py',
  'hazel/graphnet.py',
  'hazel/hsra.py',
  'hazel/io_hazel.py',
  'hazel/model.py',
  'hazel/multiprocess.py',
  'hazel/nlte_weights.py',
  'hazel/parametric.py',
  'hazel/photosphere.py',
  'hazel/rt_chromosphere.py',
  'hazel/sir.py',
  'hazel/spectrum.py',
  'hazel/stray.py',
  'hazel/tools.py',
  'hazel/transforms.py',
  'hazel/util.py',
  subdir: 'hazel',
)

py.install_sources('hazel/codes/__init__.py', subdir: 'hazel/codes')

# The neural network weights are excluded: see hazel/nlte_weights.py.
install_subdir(
  'hazel/data',
  install_dir: py.get_install_dir() / 'hazel',
  exclude_files: [
    '20211114-131045_best.prd.pth',
    '2024-09-13-11_17_34.best.pth',
  ],
)

install_data(
  'gui/hazelgui',
  install_dir: get_option('bindir'),
  install_mode: 'rwxr-xr-x',
)
