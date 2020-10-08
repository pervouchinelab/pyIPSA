import sys
from pathlib import Path
from collections import namedtuple
import pandas as pd

Record = namedtuple("Row", ("sample", "genome", "read_length", "paired", "stranded"))

out = Path(sys.argv[1])
records = []
for filename in out.rglob("*.txt"):
    sample = str(filename.parent).split("/")[-1]
    with filename.open("r") as f:
        for line in f:
            if " is " not in line:
                continue
            left, right = line.strip().split(" is ")
            if "Genome" in left:
                genome = right
            if "Read" in left:
                read_length = int(right)
            if right.endswith("-end"):
                paired = True if right[:-4] == "pair" else False
            if right.endswith("stranded"):
                stranded = True if len(right) == 8 else False
    records.append(Record(sample, genome, read_length, paired, stranded))

df = pd.DataFrame(records)
df.sort_values(by="sample").to_csv("comparison.tsv", sep="\t", index=False)
