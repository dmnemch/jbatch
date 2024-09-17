#!/bin/python3

import argparse
import subprocess 
import os.path
from pathlib import Path
import yaml

def get_parser():
    parser = argparse.ArgumentParser(
        prog='jb',
        description="Simple tool to wrap and execute sbatch commands from Jupyter Notebook cells with conda envs"
    )
    parser.add_argument("-a", "--account", metavar="<account>", help="Account to charge for resource usage")
    parser.add_argument("-r", "--reservation", metavar="<reservation>", help="Reservation name")
    parser.add_argument("-p", "--partition", metavar="<partition>", help="Partition to submit job to")
    parser.add_argument("-c", "--cpus", type=int, metavar="<cpus>", help="Number of CPUs to use")
    parser.add_argument("-m", "--mem", type=int, metavar="<memory>", help="Memory in GB (without 'G' suffix)")
    parser.add_argument("-t", "--time", type=int, metavar="<time>", help="Time in hours")
    parser.add_argument("--conda", metavar="<conda>", help="Activate conda environment by name or path")
    parser.add_argument("--logdir", metavar="<logdir>", help="Destination directory for .out and .err files")
    parser.add_argument("--name", metavar="<name>", help="Base name for log (.out & .err) files and for the sbatch job, use %cmd for first word in cmd and %j for job id")
    parser.add_argument("--config", metavar="<config>", help="Path to jb_config.yaml file")
    parser.add_argument("-v", "--verbosity", default=1, type=int, metavar="<verbosity>", help="Verbosity level: 0 - quiet, 1 - Job ID, 2 - params and cmd")
    parser.add_argument("--dry", action="store_true", help="Simulate job submission without executing commands")
    parser.add_argument("-d", "--dependency", metavar="<dependency>", help="Job dependencies")
    parser.add_argument("cmd", nargs=argparse.REMAINDER, metavar="<cmd>", help="Command to wrap into sbatch script")
    return parser


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
      print(f"\033[32mNo config file is found! Generated default config in {config_path}, \033[31mbut it lacks crucial arguments (e.g. partition) – please add them manualy \033[0m")
  if params["config"]:
    try:
      config = yaml.safe_load(Path(params["config"]).read_text())
    except:
      raise ValueError("\033[31m Config file given but loading wasn't sucsessful: check existence and format! \033[0m")
    for k,v in config.items():
      if k not in params or not params[k]: params[k] = v
  params["cmd"] = " ".join(params["cmd"])
  params["prefix"] = params["prefix"] + " && " if "prefix" in params.keys() and params["prefix"] else "" 
  params["conda"] = params["conda_prefix"] + " " + params["conda"] + " && " if "conda" in params.keys() and params["conda"] else "" 
  params["dependency"] = f'--dependency=afterok:{params["dependency"]}' if params["dependency"] else ""
  params["reservation"] = f'--reservation={params["reservation"]}' if params["reservation"] else ""
  params["account"] = f'--account {params["account"]}' if params["account"] else ""
  params["partition"] = f'--partition {params["partition"]}' if params["partition"] else ""

  return params

def run_cmd(cmd, parser):
  try:
    result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    return result
  except subprocess.CalledProcessError as e:
    parser.print_usage()
    print(f"\033[31mSbatch execution failed with error:\033[0m\n{e.stderr}")

def main():
  parser = get_parser()
  args = parser.parse_args()
  params = generate_params(vars(args))
  if not params["cmd"]:
    parser.print_usage()
    raise ValueError("\033[31m Empty command given! Void can't be submited to slurm! \033[0m")
  if params["verbosity"] > 1:
    print("\033[33mParams: \033[0m")
    for k,v in params.items():
      print("", k, v, sep="\t")
  if "%cmd" in params["name"]:
   params["name"] = params["cmd"].split()[0].join(params["name"].split("%cmd"))
  cmd2wrap = f'''sbatch {params["account"]} {params["partition"]} {params["reservation"]} {params["dependency"]} \
-o '{params["logdir"]}/{params['name']}.out' -e '{params["logdir"]}/{params["name"]}.err' --job-name={params["name"]} \
-c {params["cpus"]} --mem={params["mem"]}G --time={params["time"]}:00:00 \
--parsable --wrap "/bin/bash -c '{params["prefix"]}{params["conda"]}{params["cmd"]}'"'''
  if params["verbosity"] > 1: print("\033[33mCommand: \033[0m", cmd2wrap, sep="\n")
  if not params["dry"]:
    job_id = run_cmd(cmd2wrap, parser).stdout.strip()
    if params["verbosity"] > 0 : print(job_id)
  else:
    print("\033[33mUsing dry mode! No commands are submited to slurm!\033[0m")

if __name__ == "__main__":
    main()