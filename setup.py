from setuptools import find_packages, setup

setup(
    name='sdf',
    version='2.0.0',
    description='all your seds are belong to us',
    url='https://github.com/drgmk/sdf',
    author='Grant M. Kennedy',
    author_email='g.kennedy@warwick.ac.uk',
    license='MIT',
    packages=find_packages(exclude=('tests', 'tests.*')),
    classifiers=['Programming Language :: Python :: 3'],
    install_requires=[
        'astropy', 'bokeh', 'corner',
        'extinction', 'emcee', 'filelock', 'bleach',
        'jinja2', 'matplotlib', 'mysql-connector-python', 'numpy',
        'pymultinest', 'pytest', 'PyYAML', 'requests', 'scipy'
        ],
    extras_require={
        'docs': ['myst-parser', 'numpydoc', 'sphinx'],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': ['sdf-fit=sdf.scripts:sdf_fit',
                            'sdf-fits=sdf.scripts:sdf_fits',
                            'sdf-sample=sdf.scripts:sdf_sample',
                            'sdf-index=sdf.scripts:sdf_index',
                            'sdf-jobs=sdf.scripts:sdf_jobs',
                            'sdf-cleanup=sdf.scripts:sdf_cleanup']
        },
    zip_safe=False
    )
