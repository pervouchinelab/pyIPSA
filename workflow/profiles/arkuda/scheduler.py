#!/usr/bin/env python3
import os, sys
import logging, traceback

logging.basicConfig(
    level=logging.INFO,
    format="CLUSTER: %(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

## Beginning of the script

from subprocess import Popen, PIPE
import yaml


os.makedirs("cluster_log", exist_ok=True)


# let snakemake read job_properties
from snakemake.utils import read_job_properties

wrapper_directory = os.path.dirname(__file__)

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# default paramters defined in cluster_spec (accessed via snakemake read_job_properties)
cluster_param = job_properties["cluster"]

if job_properties["type"] == "single":
    cluster_param["name"] = job_properties["rule"]
elif job_properties["type"] == "group":
    cluster_param["name"] = job_properties["groupid"]
else:
    raise NotImplementedError(
        f"Don't know what to do with job_properties['type']=={job_properties['type']}"
    )


# don't overwrite default parameters if defined in rule (or config file)
if ("threads" in job_properties) and ("threads" not in cluster_param):
    cluster_param["threads"] = job_properties["threads"]

for res in ["time_min", "mem_mb"]:
    if (res in job_properties["resources"]) and (res not in cluster_param):
        cluster_param[res] = job_properties["resources"][res]


# check which system you are on and load command command_options
key_mapping_file = os.path.join(wrapper_directory, "key_mapping.yaml")
command_options = yaml.load(open(key_mapping_file), Loader=yaml.BaseLoader)
system = command_options["system"]
command = command_options[system]["command"]

key_mapping = command_options[system]["key_mapping"]

# construct command:
for key in cluster_param:
    if key not in key_mapping:
        logging.warning(
            f"parameter '{key}' not in keymapping! It would be better if you add the key to the file: {key_mapping_file} \n I try without the key!"
        )
    else:
        command += " "
        command += key_mapping[key].format(cluster_param[key])

command += " {}".format(jobscript)

logging.info("submit command: " + command)

p = Popen(command.split(), stdout=PIPE, stderr=PIPE)
output, error = p.communicate()
if p.returncode != 0:

    error_message = (
        "Job can't be submitted\n" + output.decode("utf-8") + error.decode("utf-8")
    )

    logging.error(error_message)

    # give it to snakemake
    print(error_message)

    exit(p.returncode)
else:
    res = output.decode("utf-8")

    if system == "lsf":
        import re

        match = re.search(r"Job <(\d+)> is submitted", res)
        jobid = match.group(1)

    elif system == "pbs" or system == "pbspro":
        jobid = res.strip().split(".")[0]

    else:
        jobid = int(res.strip().split()[-1])

    print(jobid)
