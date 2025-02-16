from distutils.core import setup, Extension

all_macros = [('MOTION_DEBUG', None), ('PYTHON', None)]


rangetree = Extension('ring_rangetree',
                    sources = ['c/ring_rangetree.c'],
                    extra_compile_args = ["-O3"],
                    define_macros=all_macros)

setup (name = 'ring_rangetree',
       version = '1.0',
       description = 'rangetree with loopy dimensions',
       ext_modules = [rangetree])
