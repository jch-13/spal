from setuptools import setup, find_packages

setup(
    name='spal',
    version='0.1.0',
    packages=find_packages(),
    url='',
    license='',
    author='Cesare de Filippo',
    author_email='casare_filippo@eva.mpg.de',
    description='Report the estimates of spurious alignments as well as the the count of spurious and true alignments',
    install_requires=['pysam', 'scipy'],
    python_requires='>=3.5',
    entry_points={'console_scripts': ['spal=spal:main']},
)
