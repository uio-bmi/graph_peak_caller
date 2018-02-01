from setuptools import setup

setup(name='graph_peak_caller',
      version='1.0.0',
      description='Offset based graph',
      url='http://github.com/uio-bmi/graph_peak_caller',
      author='Ivar Grytten and Knut Rand',
      author_email='',
      license='MIT',
      packages=['graph_peak_caller'],
      zip_safe=False,
      install_requires=['pymysql', 'numpy', 'filecache', 'scipy', 'pybedtools', 'pyBigWig',
                        'memory_profiler', 'python-coveralls', 'matplotlib', 'codecov', 'pytest-cov'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ]

      )
