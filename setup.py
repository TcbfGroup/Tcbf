from setuptools import setup,Extension

setup(name = "tcbf",
      description='Tcbf identify the conservative TAD boundary between multipy species',
      author = "He Xin",
      license= 'MIT License',
      license_files = ('LICENSE',),
      packages = ["tcbf"],
      url = "https://github.com/hexin010101/Tcbf",
     scripts=['bin/tcbf','bin/tcbf_syn_process'])