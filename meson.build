project('lattice-cpp', 'cpp', default_options : ['cpp_std=c++2a'])

source = files(['src/main.cpp',
                'src/parse.cpp',
                'src/search.cpp',
                'src/matrix.cpp',
                'src/optimize.cpp',
                'src/format.cc'])

threads = dependency('threads')
gmp = dependency('gmp')

exe = executable('lattice', source, dependencies: [threads, gmp], install: true, install_dir: join_paths(meson.current_source_dir(), 'out/'))
