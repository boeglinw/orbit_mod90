project('fluxpy',
        ['c', 'fortran'],
        version : '0.1',
        meson_version: '>= 1.1.0',
        default_options : [
                            'warning_level=1',
                            'buildtype=release'
                          ])
fc = meson.get_compiler('fortran')

py = import('python').find_installation('''/opt/anaconda3/envs/cern_root/bin/python''', pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

inc_np = include_directories(incdir_numpy, incdir_f2py)
# gh-25000
quadmath_dep = fc.find_library('quadmath', required: false)

flux = declare_dependency(link_args : ['-lflux'])
read_eqdsk = declare_dependency(link_args : ['-lread_eqdsk'])
modules = declare_dependency(link_args : ['-lmodules'])
spline = declare_dependency(link_args : ['-lspline'])
lib_dir_0 = declare_dependency(link_args : ['''-L/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90//flux/'''])
lib_dir_1 = declare_dependency(link_args : ['''-L/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90//read_eqdsk/'''])
lib_dir_2 = declare_dependency(link_args : ['''-L/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90//include/'''])
lib_dir_3 = declare_dependency(link_args : ['''-L/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90//spline/'''])

py.extension_module('fluxpy',
                     [
                     '''fluxpy.f90''',
                     '''fluxpymodule.c''',
                     '''fluxpy-f2pywrappers2.f90''',
                     fortranobject_c
                     ],
                     include_directories: [
                     inc_np,

                     ],
                     dependencies : [
                     py_dep,
                     quadmath_dep,

                     flux,
                                          read_eqdsk,
                                          modules,
                                          spline,
                     lib_dir_0,
                                          lib_dir_1,
                                          lib_dir_2,
                                          lib_dir_3,
                     ],

                     install : true)
