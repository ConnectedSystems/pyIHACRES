try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup
    import re
    import os

    def find_packages(path='.'):
        ret = []
        for root, dirs, files in os.walk(path):
            if '__init__.py' in files:
                ret.append(re.sub('^[^A-z0-9_]+', '', root.replace('/', '.')))
        return ret
    # End find_packages()
# End try

pkg_list = ['numpy', 'pandas', 'matplotlib', 'seaborn', 'PyYAML', 'sklearn']
config = {
    'description': 'pyIHACRES',
    'author': 'Takuya Iwanaga, Barry F.W. Croke',
    'url': 'https://github.com/ConnectedSystems/pyIHACRES',
    'download_url': 'https://github.com/ConnectedSystems/pyIHACRES/archive/master.zip',
    'author_email': 'iwanaga.takuya@anu.edu.au',
    'version': '0.1',
    'install_requires': pkg_list,
    'packages': find_packages(),
    'scripts': [],
    'name': 'pyIHACRES'
}

setup(**config)