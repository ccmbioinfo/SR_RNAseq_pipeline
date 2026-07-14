import os
import pandas as pd
import math
import re

tissue=config["sampleInfo"]["tissue"]
sampleID=config["sampleInfo"]["sampleID"]
output_dir=config["sampleInfo"]["dir"]
sample_config=os.path.join(config["sampleInfo"]["dir"], "config.yaml")

trimming=True if config["script_parameters"]["trimming"] == 'TRUE' else False
