from setuptools import setup, find_packages

# versioning
import versioneer

setup(
    name='topobuilder',
    version=versioneer.get_version(),

    description='generation of tailored protein topologies',
    long_description='topobuilder is a python library that, in combination with the '
                     'Rosetta suite, aids in the semi-automatic generation of user-defined '
                     'topologies. It expressely adds the capability to tailor those topologies '
                     'around a functional motif of interest.',

    # The project's main homepage.
    url='https://github.com/jaumebonet/topobuilder',

    # Author details
    author='Jaume Bonet',
    author_email='jaume.bonet@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],

    project_urls={
        'Documentation': 'http://jaumebonet.cat/topobuilder',
        'Source': 'https://github.com/jaumebonet/topobuilder/',
        'Tracker': 'https://github.com/jaumebonet/topobuilder/issues',
    },

    platforms='UNIX',
    keywords='development',

    install_requires=['numpy', 'scipy', 'networkx', 'svgwrite', 'transforms3d', 'bottle'],

    packages=find_packages(exclude=['docs', 'demo', 'sphinx-docs']),
    include_package_data=True,
    package_data={
        'topobuilder': []
    },
    entry_points={
        'console_scripts':
            ['topo.case=topobuilder.io.case:case_template',
             ]
    },
    cmdclass=versioneer.get_cmdclass(),
)
