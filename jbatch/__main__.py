#!/bin/python3

import argparse
import subprocess 
import os.path
from pathlib import Path
import yaml

def parse_args():
  parser = argparse.ArgumentParser(
    prog='jb',
    description="\033[35;1m Simple tool to wrap and execute sbatch commands from Jupyter Notebook cells with conda envs \033[0m")
  parser.add_argument("-a", "--account", metavar="", help="Account to charge for resource usage :)")
  parser.add_argument("-r", "--reservation", metavar="", help="Reservation name")
  parser.add_argument("-p", "--partition", metavar="", help="Partition to submit job to")
  parser.add_argument("-c", "--cpus", type=int, metavar="", help="Number of CPUs to use")
  parser.add_argument("-m", "--mem", type=int, metavar="", help="Memory in GB (without 'G' suffix)")
  parser.add_argument("-t", "--time", type=int, metavar="", help="Time in hours")
  parser.add_argument("--conda", metavar="", help="Activate conda environment by name or path")
  parser.add_argument("--logdir", metavar="", help="Destination directory for .out and .err files")
  parser.add_argument("--name", metavar="", help="Base name for .out and .err files")
  parser.add_argument("--config", metavar="", help="Path to jb_config.yaml file")
  parser.add_argument("-v", "--verbosity", type=int, metavar="", help="Verbosity level: 0 - quiet, 1 - Job ID, 2 - params and cmd")
  parser.add_argument("--dry", action="store_true", help="Simulate job submission without executing commands")
  parser.add_argument("-d", "--dependency", metavar="", help="Job dependencies")
  parser.add_argument("cmd", nargs=argparse.REMAINDER, metavar="", help="Command to wrap into sbatch script")
  return parser.parse_args()

def generate_config(config_path):
  default_config = {'account': '', 'partition': '', 'reservation': '', 'cpus': 4, 'mem': 4, 'time': 24, 'prefix': '', 'conda_prefix': 'conda activate', 'logdir': '.', 'name': '%j', 'verbosity': 1}
  config_dir = os.path.dirname(config_path)
  if not os.path.exists(config_dir):
    os.makedirs(config_dir)
  with open(config_path, 'w') as out:
    yaml.dump(default_config, out, default_flow_style=False)

def generate_params(params):
  """Takes cli arguments and config file and generates consensus params"""
  if not params["config"]:
    if os.path.exists("jb_config.yaml"):
      if params["verbosity"] > 1 : print()
      params["config"] = "jb_config.yaml"
    elif os.path.exists(f"{os.path.expanduser('~')}/.config/jb/jb_config.yaml"): 
      params["config"] = f"{os.path.expanduser('~')}/.config/jb/jb_config.yaml"
    else:
      config_path = f"{os.path.expanduser('~')}/.config/jb/jb_config.yaml"
      generate_config(config_path)
      print(f"\033[32mNo config file is found! Generated default config in {config_path}, \033[31mbut it lacks crucial argument (e.g. partition) – please add them manualy \033[0m")
  if params["config"]:
    try:
      config = yaml.safe_load(Path(params["config"]).read_text())
    except:
      raise ValueError("Config file given but loading wasn't sucsessful: check format!")
    for k,v in config.items():
      if k not in params or not params[k]: params[k] = v
  params["cmd"] = " ".join(params["cmd"])
  params["prefix"] = params["prefix"] + " && " if params["prefix"] else "" 
  params["conda"] = params["conda_prefix"] + " " + params["conda"] + " && " if params["conda"] else "" 
  params["dependency"] = f'--dependency=afterok:{params["dependency"]}' if params["dependency"] else ""
  params["reservation"] = f'--reservation={params["reservation"]}' if params["reservation"] else ""
  return params

def run_cmd(cmd):
  try:
    result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    return result
  except subprocess.CalledProcessError as e:
    print(f"Sbatch execution failed with error:\n{e.stderr}")

def main():
  args = parse_args()
  params = generate_params(vars(args))
  if not params["cmd"]:
    raise ValueError("\033[31;1m No command are submited to slurm!\033[0m")
  if params["verbosity"] > 1:
    print("\033[33;1mParams: \033[0m")
    for k,v in params.items():
      print("", k, v, sep="\t")
  cmd = f'''sbatch --account {params["account"]} --partition {params["partition"]} \
{params["reservation"]} {params["dependency"]} \
-o '{params["logdir"]}/{params['name']}.out' -e '{params["logdir"]}/{params["name"]}.err' \
-c {params["cpus"]} --mem={params["mem"]}G --time={params["time"]}:00:00 \
--parsable --wrap "/bin/bash -c '{params["prefix"]}{params["conda"]}{params["cmd"]}'"'''
  if params["verbosity"] > 1: print("\033[33;1mCommand: \033[0m", cmd, sep="\n")
  if not params["dry"]:
    job_id = run_cmd(cmd).stdout.strip()
    if params["verbosity"] > 0 : print(job_id)
  else:
    print("\033[33mUsing dry mode! No commands are submited to slurm!\033[0m")
    if params["verbosity"] < 2: print("To print params and sbatch command, set --verbosity to 2")

if __name__ == "__main__":
    main()