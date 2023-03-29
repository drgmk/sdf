from setuptools import setup

setup(
    name='sdf',
    version='0.1',
    description='all your seds are belong to us',
    url='https://github.com/drgmk/sdf',
    author='Grant M. Kennedy',
    author_email='g.kennedy@warwick.ac.uk',
    license='MIT',
    packages=['sdf'],
    classifiers=['Programming Language :: Python :: 3'],
    install_requires=[
        'astropy', 'bokeh', 'corner',
        'extinction', 'emcee', 'filelock', 'bleach',
        'jinja2', 'matplotlib', 'mysql-connector-python', 'numpy',
        'pymultinest', 'pytest', 'requests', 'scipy'
        ],
    include_package_data=True,
    entry_points={
        'console_scripts': ['sdf-fit=sdf.scripts:sdf_fit',
                            'sdf-sample=sdf.scripts:sdf_sample',
                            'sdf-cleanup=sdf.scripts:sdf_cleanup']
        },
    zip_safe=False
    )
