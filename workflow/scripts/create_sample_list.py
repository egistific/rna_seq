import os
import yaml

from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[2]
print(PROJECT_DIR)

with open(os.path.join(PROJECT_DIR, "config/config.yaml")) as config_file:
    config = yaml.safe_load(config_file)

READS_DIR = config["data"]["reads_dir"]
METADATA_FILEPATH = os.path.join(PROJECT_DIR, config["data"]["metadata"])
SAMPLES_FILEPATH = os.path.join(PROJECT_DIR, config["data"]["samples"])

with open(SAMPLES_FILEPATH, "w") as samples:
    samples.write("sample_id\tdonor\tcondition\tfq1\tfq2\n")    

fastq_filepaths = list(Path(os.path.join(PROJECT_DIR, READS_DIR)).glob("*.fastq.gz"))

with open(METADATA_FILEPATH, "rt") as metadata, open(SAMPLES_FILEPATH, "a") as samples:
    for line in metadata:
        text = line.strip().split(",")
        gsm, sample_id, run = text
        # use only "untreated" and "Dex" groups for now
        if "Alb" in sample_id:
            continue
        donor, condition = sample_id.split("_", maxsplit=1)
        fq1 = [f"{READS_DIR}/{f.name}" for f in fastq_filepaths if run in f.name and "_1" in f.name]
        fq2 = [f"{READS_DIR}/{f.name}" for f in fastq_filepaths if run in f.name and "_2" in f.name]
        assert len(fq1) == 1 and len(fq2) == 1, "Check fastq files!"
        samples.write(f"{sample_id}\t{donor}\t{condition}\t{fq1[0]}\t{fq2[0]}\n")

print("Created sample list!")


