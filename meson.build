project('lattice-cpp', 'cpp', default_options : ['cpp_std=c++2a'])

inc = include_directories('include')
src = files(['src/la.cpp'])
exe = executable('lattice', src, include_directories: inc, install: true, install_dir: join_paths(meson.current_source_dir(), 'out/'))
