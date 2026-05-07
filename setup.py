from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        'ODdoubleR',
        sources=[
            'src/ODdoubleR/bindings.cpp',
            'src/ODdoubleR/main.cpp',
            'src/ODdoubleR/init_det.cpp'
        ],
        include_dirs=[
            pybind11.get_include(),
            pybind11.get_include(user=True),
            'include/'
        ],
        language='c++',
        extra_compile_args=['-std=c++17', '-O3'],
    ),
]

setup(
    name='ODdoubleR',
    version='1.0.0',
    author='Alex_kiam',
    author_email='sasha.gafarov.03@mail.ru',
    description='Double R determination module',
    long_description_content_type='text/markdown',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.10'],
    python_requires='>=3.7',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
    ],
    zip_safe=False,
)