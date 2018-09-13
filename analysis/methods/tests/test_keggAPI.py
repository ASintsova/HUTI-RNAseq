import sys
sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/')
from methods import keggAPI

from Bio.KEGG import Enzyme
from Bio.KEGG import REST



def _write_kegg(item, info, indent=12):

      s = ""
      for line in info:
          partial_lines = line.splitlines()
          for l in partial_lines:
              s = s + item.ljust(indent) + l + "\n"
              if item is not "":  # ensure item is only written on first line
                  item = ""
      return s


print(_write_kegg("ENTRY", ["b0356"]))