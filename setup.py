from setuptools import setup

setup(name = "tcbf",
      description='Tcbf used to identify the conservative TAD boundary between multipy species',
      author = "He Xin",
      license= 'MIT License',
      license_files = ('LICENSE',),
      packages = ["tcbf"],
      include_package_data=True,
      url = "https://github.com/hexin010101/Tcbf",
     scripts=['bin/tcbf','bin/tcbf_syn_process','bin/plot_TAD_bound_synteny'])