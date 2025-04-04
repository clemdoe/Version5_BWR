#
# python3 setup_cle2000.py install --home=.
#
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib

def main():
  import os
  from distutils.sysconfig import get_config_var
  incdir = os.path.join(get_python_lib(plat_specific=1), "numpy/core/include")
  mach = os.path.basename(os.getcwd())
  Code = os.environ.get("CODE_EMBEDDED", None) # Code selection
  FortranLib = os.environ.get("FORTRANPATH", None) # directory with libgfortran.a
  HDF5Lib = os.environ.get("HDF5_API", None) # directory with libhdf5.a
  pylib = os.path.basename(get_config_var("LIBDIR")) # get lib or lib64
  print("install Cle2000 binding to", Code, "on directory",mach, "pylib=",pylib)

  if Code == "GANLIB":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with GANLIB",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach,FortranLib,HDF5Lib],
                     runtime_library_dirs=[HDF5Lib],
                     libraries=["Ganlib","gfortran","hdf5"] ) ])
  elif Code == "TRIVAC":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with TRIVAC",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__trivac__', None)],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach,FortranLib,HDF5Lib,
                     "../../../Utilib/lib/"+mach,"../../../Trivac/lib/"+mach],
                     runtime_library_dirs=[HDF5Lib],
                     libraries=["Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "DRAGON":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with DRAGON",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__dragon__', None)],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach,FortranLib,HDF5Lib,"../../../Utilib/lib/"+mach
                     ,"../../../Trivac/lib/"+mach,"../../../Dragon/lib/"+mach],
                     runtime_library_dirs=[HDF5Lib],
                     libraries=["Dragon","Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "DONJON":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with DONJON",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__donjon__', None)],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach,FortranLib,HDF5Lib,
                     "../../../Utilib/lib/"+mach,"../../../Trivac/lib/"+mach,
                     "../../../Dragon/lib/"+mach,"../../../Donjon/lib/"+mach],
                     runtime_library_dirs=[HDF5Lib],
                     libraries=["Donjon","Dragon","Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "GANLIB_OMP":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with GANLIB_OMP",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     extra_f90_compile_args = ["-fopenmp"],
                     extra_link_args = ["-lgomp"],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach,FortranLib,HDF5Lib],
                     runtime_library_dirs=[HDF5Lib],
                     libraries=["Ganlib","gfortran","hdf5"] ) ])
  elif Code == "TRIVAC_OMP":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with TRIVAC_OMP",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__trivac__', None)],
                     extra_f90_compile_args = ["-fopenmp"],
                     extra_link_args = ["-lgomp"],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach,FortranLib,HDF5Lib,
                     "../../../Utilib/lib/"+mach,"../../../Trivac/lib/"+mach],
                     runtime_library_dirs=[HDF5Lib],
                     libraries=["Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "DRAGON_OMP":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with DRAGON_OMP",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__dragon__', None)],
                     extra_f90_compile_args = ["-fopenmp"],
                     extra_link_args = ["-lgomp"],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach,FortranLib,HDF5Lib,"../../../Utilib/lib/"+mach,
                     "../../../Trivac/lib/"+mach,"../../../Dragon/lib/"+mach],
                     runtime_library_dirs=[HDF5Lib],
                     libraries=["Dragon","Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "DONJON_OMP":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with DONJON_OMP",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__donjon__', None)],
                     extra_f90_compile_args = ["-fopenmp"],
                     extra_link_args = ["-lgomp"],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach,FortranLib,HDF5Lib,
                     "../../../Utilib/lib/"+mach,"../../../Trivac/lib/"+mach,
                     "../../../Dragon/lib/"+mach,"../../../Donjon/lib/"+mach],
                     runtime_library_dirs=[HDF5Lib],
                     libraries=["Donjon","Dragon","Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  else:
    raise ValueError(Code+" is not implemented for distutils bindings")
if __name__ == "__main__":
  main()
