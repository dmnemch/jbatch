import subprocess
import sys

args = sys.argv[1:]

def run_bash_script(script, args):
  try:
    result = subprocess.run(['bash', script] + args, check=True, text=True, capture_output=True)
    print(result.stdout)
  except subprocess.CalledProcessError as e:
    print(f"An error occurred while running the script: {e}")
    print(e.stdout, e.stderr)

def err():
  run_bash_script("./err", args)

def out():
  run_bash_script("./out", args)