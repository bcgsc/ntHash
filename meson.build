project(
    'ntHash',
    'cpp',
    version: '2.4.0',
    license: 'MIT',
    default_options: ['cpp_std=c++17', 'werror=true', 'warning_level=3'],
)

include_dirs = [include_directories('include')]

sources = ['src/kmer.cpp', 'src/seed.cpp']

nthash_lib = static_library(
    'nthash',
    sources,
    include_directories: include_dirs,
    install: true,
    install_dir: 'lib',
)

lib_dep = declare_dependency(
    link_with: nthash_lib,
    include_directories: include_dirs,
)

if get_option('buildtype') != 'release'
    executable(
        'bench',
        'examples/benchmark.cpp',
        dependencies: [lib_dep],
        include_directories: include_dirs,
    )
    test_exec = executable(
        'nthash',
        'tests/tests.cpp',
        dependencies: [lib_dep],
        include_directories: include_dirs,
        install: true,
    )
    test('tests', test_exec)
endif
