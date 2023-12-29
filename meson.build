project('SHGOFAST','cpp','cython', version : '1.0.0',
        default_options : ['optimization=3',
                           'cpp_std=c++20','build.cpp_std=c++20'])

add_project_arguments(['-march=native'], language : 'cpp')

include = include_directories('include')

python = import('python').find_installation()
dep_py = python.dependency()

incpath_np = run_command(
  'gravelspoon',
  check : true
).stdout().strip()
incdir_np = include_directories(incpath_np)

python.extension_module('vertex',
                        'shgofast/vertex.pyx',
                        dependencies : [dep_py],
                        include_directories : [include, incdir_np],
                        override_options : ['cython_language=cpp'])