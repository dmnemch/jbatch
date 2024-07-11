## Jupyter Notebooks 🤝 Conda 🤝 Slurm
An alternative for writing separate sbatch scripts and bash pipelines for them – `jb` is a simple tool to wrap and execute sbatch comands from Juputer Notebook cells with conda envs:
 - ⚙️ Config files with default sbatch parametres – e.g. `--account`, `--partition`
 - 🚀 Python variables inside sbatch commands – e.g. `for sample in samples: ...`, see example №1
 - 🐍 Conda envs with `--conda` flag, see example №2 
 - 🔗 Manage dependences for jobs, see example №4

## Arguments
| Argument                | Description                                                                                                 |
| ----------------------- | ----------------------------------------------------------------------------------------------------------- |
| `--account`, `-a`        | Account to charge for resource usage                                                                        |
| `--reservation`, `-r`   | Reservation name                                                                                            |
| `--partition`, `-p`     | Partition to submit the job to                                                                              |
| `--cpus`, `-c`          | Number of CPUs to use (e.g. 4)                                                                              |
| `--mem`, `-m`           | Memory in GB without 'G' suffix (e.g. 8)                                                                    |
| `--time`, `-t`          | Time in hours (e.g. 24)                                                                                     |
| `--conda`               | Activate conda environment by name or path. Conda executable must be in the $PATH                           |
| `--logdir`              | Destination directory for .out and .err files                                                               |
| `--name`                | Base name for .out and .err files                                                                           |
| `--config`              | Path to jb_config.yaml file. If no config given, trying to find one in the current dir and in ~/.config/jb  |
| `--verbosity`, `-v`     | Verbosity level: 0 - quiet, 1 - Job ID, 2 - params and cmd                                                  |
| `--dry`                 | Simulate job submission without executing commands                                                          |
| `--dependency`, `-d`    | Job dependencies                                                                                            |

## Config
Config data in `yaml` fromat may be located in:
  - ~/.config/jb/jb_config.yaml [default]
  - ./jb_config.yaml
  - anywhere if provided with `--config` flag

Any argument can be specified here and there are two extra options:
`prefix` for a command to be executed before anything else, including conda actiavtion

`conda_prefix` for a command to be executed before any conda env activation (some clusters require loading anaconda module first)

*Note*: no '&&' are needed in this prefixes – they will be added by `jb`

## Examples
#### 📌 Example №0: Simple execution of sbatch command
```python
! jb bwa index ~/ref/hg38.fa
```

#### 🚀 Example №1: Using python variables
```python
import glob
ref = "~/ref/hg38.fa" 
for sample in glob.glob("~/data/*.fq"):
  ! jb bwa mem {ref} {sample} -o {sample[:-2] + "sam"}
```

#### 🐍 Example №2: Using conda environment
```python
path = "~/data/sample1.sam"
! jb --conda ngs samtools view {path}
```

#### ⚠️ Example №3: Using pipes and redirections
*Tip*: To use special bash symbols like "|", ">" or "&" use escape symbol "\"
```python
! jb echo "hi!" \> test123 \&\& sleep 5
```

#### 🔗 Example №4: Managing dependencies
*Tip*: Use `--dry` and `-v 2` to show command to execute without execution (but dependencies aren't shown)
```python
job_ids = []
job_id = ! jb echo "Hello" \> "hi.txt"
job_ids.append(job_id)
job_id = ! jb echo "World" \>\> "hi.txt"
job_ids.append(job_id)

job_ids = ",".join([_ for __ in job_ids for _ in __])

! jb -d {job_ids} cat "hi.txt"
```

## Future features (maybe...)
 - 🪵 Profile executes comands for resource usage (cpus, memory, disc i/o, time)
 - 🏛️ History of used commands
 - ✨ Fancy dependencies
 - ➕ Extra arguments for sbatch
