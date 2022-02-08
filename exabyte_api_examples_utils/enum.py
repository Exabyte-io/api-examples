# TODO: better name? enum used elsewhere in exabyte codebase but not sure why

from os.path import abspath, dirname, join

PORT = 443
SECURE = True
VERSION = "2018-10-01"
HOST = "platform.exabyte.io"

module_dir = dirname(abspath(__file__))

gamma_alumina_file = join(module_dir, "assets/gamma_alumina_digne_et_al.poscar")
mp_978534_file = join(module_dir, "assets/mp-978534.poscar")
