from setuptools import setup

setup(
    name='sdf',
    version='0.1',
    description='all your seds are belong to us',
    url='http://github.com/drgmk/sdf',
    author='Grant M. Kennedy',
    author_email='g.kennedy@warwick.ac.uk',
    license='MIT',
    packages=['sdf'],
    classifiers=['Programming Language :: Python :: 3'],
    install_requires = [
        'astropy >= v2.0.0','binarytree','bokeh','corner',
        'extinction','emcee','filelock',
        'jinja2','matplotlib','mysql-connector','numpy',
        'pymultinest','requests','scipy'
        ],
    entry_points = {
        'console_scripts': ['sdf-fit=sdf.scripts:sdf_fit',
                            'sdf-sample=sdf.scripts:sdf_sample',
                            'sdf-cleanup=sdf.scripts:sdf_cleanup']
        },
    zip_safe=False
    )
