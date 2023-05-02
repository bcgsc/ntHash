project(
    'ntHash',
    'cpp',
    version: '2.4.0',
    license: 'MIT',
    default_options: ['cpp_std=c++17', 'werror=true', 'warning_level=3'],
)

include_dirs = [
    include_directories('include')
]

sources = ['src/nthash.cpp', 'src/lowlevel.cpp']

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
    catch2 = dependency(
        'Catch2',
        method: 'cmake',
        modules: ['Catch2::Catch2WithMain'],
        required: false,
    )
    test_exec = executable(
        'nthash',
        'tests/tests.cpp',
        dependencies: [lib_dep, catch2],
        include_directories: include_dirs,
        install: true,
    )
    test('tests', test_exec)
endif
