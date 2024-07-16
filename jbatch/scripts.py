import os
import sys
from importlib.resources import files

def run_bash_script(script):
  resource = files("jbatch")
  os.environ["PATH"] = f"{resource}:{os.getenv('PATH')}"
  script = files("jbatch").joinpath(script)
  print(script)
  os.execv(script, sys.argv)

def err():
  run_bash_script("err")

def out():
  run_bash_script("out")

if __name__ == "__main__":
  err()