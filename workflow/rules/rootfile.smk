import os
from glob import glob
import re

# Directories
BUILD_DIR = "build"
EEEMCal_DATA_DIR = "data/DESY_2025/data/beam"
FoCal_DATA_DIR = "data/SPS_2024/data/beam"
GENERAL_DATA_DIR = "data"
DUMP_DIR = "dump"
LOG_DIR = "logs"

# Find all RunXXX.h2g files
data_files_eeemcal = glob(os.path.join(EEEMCal_DATA_DIR, "Run*.h2g"))
data_files_focal = glob(os.path.join(FoCal_DATA_DIR, "Run*.h2g"))

# Extract RunXXX identifiers
run_ids_eeemcal = [re.search(r'Run(\d+)\.h2g$', f).group(1) for f in data_files_eeemcal if re.search(r'Run(\d+)\.h2g$', f)]
run_ids_focal = [re.search(r'Run(\d+)\.h2g$', f).group(1) for f in data_files_focal if re.search(r'Run(\d+)\.h2g$', f)]

# Rule to process data files using the compiled 100_Rootifier script
# rule Rootifier_EEEMCal:
#     input:
#         BUILD_DIR + "/100_Rootifier",
#         EEEMCal_DATA_DIR + "/Run{run_id}.h2g"
#     output:
#         temp(DUMP_DIR + "/100_Rootifier/EEEMCal/Run{run_id}_rootified.root")
#     log:
#         LOG_DIR + "/100_Rootifier/EEEMCal/Run{run_id}.log"
#     shell:
#         """
#         mkdir -p {LOG_DIR}/100_Rootifier/EEEMCal
#         mkdir -p {DUMP_DIR}/100_Rootifier/EEEMCal
#         {input[0]} -f {input[1]} -o {output} &> {log}
#         """

# rule Rootifier_FoCal:
#     input:
#         BUILD_DIR + "/100_Rootifier",
#         FoCal_DATA_DIR + "/Run{run_id}.h2g"
#     output:
#         temp(DUMP_DIR + "/100_Rootifier/FoCal/Run{run_id}_rootified.root")
#     log:
#         LOG_DIR + "/100_Rootifier/FoCal/Run{run_id}.log"
#     shell:
#         """
#         mkdir -p {LOG_DIR}/100_Rootifier/FoCal
#         mkdir -p {DUMP_DIR}/100_Rootifier/FoCal
#         {input[0]} -f {input[1]} -o {output} &> {log}
#         """

# rule EventRecon_EEEMCal:
#     input:
#         BUILD_DIR + "/101_EventRecon",
#         DUMP_DIR + "/100_Rootifier/EEEMCal/Run{run_id}_rootified.root"
#     output:
#         temp(DUMP_DIR + "/101_EventRecon/EEEMCal/Run{run_id}_recon.root")
#     log:
#         LOG_DIR + "/101_EventRecon/EEEMCal/Run{run_id}.log"
#     shell:
#         """
#         mkdir -p {LOG_DIR}/101_EventRecon/EEEMCal
#         mkdir -p {DUMP_DIR}/101_EventRecon/EEEMCal
#         {input[0]} -f {input[1]} -o {output} &> {log}
#         """

# rule EventRecon_FoCal:
#     input:
#         BUILD_DIR + "/101_EventRecon",
#         DUMP_DIR + "/100_Rootifier/FoCal/Run{run_id}_rootified.root"
#     output:
#         temp(DUMP_DIR + "/101_EventRecon/FoCal/Run{run_id}_recon.root")
#     log:
#         LOG_DIR + "/101_EventRecon/FoCal/Run{run_id}.log"
#     shell:
#         """
#         mkdir -p {LOG_DIR}/101_EventRecon/FoCal
#         mkdir -p {DUMP_DIR}/101_EventRecon/FoCal
#         {input[0]} -f {input[1]} -o {output} &> {log}
#         """

# rule EventMatch_EEEMCal:
#     input:
#         BUILD_DIR + "/102_EventMatch",
#         DUMP_DIR + "/101_EventRecon/EEEMCal/Run{run_id}_recon.root"
#     output:
#         DUMP_DIR + "/102_EventMatch/EEEMCal/Run{run_id}_matched.root"
#     log:
#         LOG_DIR + "/102_EventMatch/EEEMCal/Run{run_id}.log"
#     shell:
#         """
#         mkdir -p {LOG_DIR}/102_EventMatch/EEEMCal
#         mkdir -p {DUMP_DIR}/102_EventMatch/EEEMCal
#         {input[0]} -f {input[1]} -o {output} &> {log}
#         """

# rule EventMatch_FoCal:
#     input:
#         BUILD_DIR + "/102_EventMatch",
#         DUMP_DIR + "/101_EventRecon/FoCal/Run{run_id}_recon.root"
#     output:
#         DUMP_DIR + "/102_EventMatch/FoCal/Run{run_id}_matched.root"
#     log:
#         LOG_DIR + "/102_EventMatch/FoCal/Run{run_id}.log"
#     shell:
#         """
#         mkdir -p {LOG_DIR}/102_EventMatch/FoCal
#         mkdir -p {DUMP_DIR}/102_EventMatch/FoCal
#         {input[0]} -f {input[1]} -o {output} &> {log}
#         """

# rule rootify_eeemcal_all:
#     input:
#         expand(DUMP_DIR + "/102_EventMatch/EEEMCal/Run{run_id}_matched.root", run_id=run_ids_eeemcal)
#     output:
#         DUMP_DIR + "/all_matched_eeemcal.done"
#     shell:
#         """
#         touch {output}
#         """

# rule rootify_focal_all:
#     input:
#         expand(DUMP_DIR + "/102_EventMatch/FoCal/Run{run_id}_matched.root", run_id=run_ids_focal)
#     output:
#         DUMP_DIR + "/all_matched_focal.done"
#     shell:
#         """
#         touch {output}
#         """

rule Rootifier_10G:
    input:
        BUILD_DIR + "/103_Rootifier_10G",
        GENERAL_DATA_DIR + "/Run{run_id}.h2g"
    output:
        DUMP_DIR + "/103_Rootifier_10G/Run{run_id}_rootified.root"
    log:
        LOG_DIR + "/103_Rootifier_10G/Run{run_id}.log"
    shell:
        """
        mkdir -p {LOG_DIR}/103_Rootifier_10G
        mkdir -p {DUMP_DIR}/103_Rootifier_10G
        {input[0]} -f {input[1]} -o {output} &> {log}
        """

rule EventRecon10G:
    input:
        BUILD_DIR + "/101_EventRecon",
        DUMP_DIR + "/103_Rootifier_10G/Run{run_id}_rootified.root"
    output:
        DUMP_DIR + "/101_EventRecon/Run{run_id}_recon.root"
    log:
        LOG_DIR + "/101_EventRecon/Run{run_id}_recon.log"
    shell:
        """
        mkdir -p {LOG_DIR}
        mkdir -p {DUMP_DIR}/101_EventRecon
        {input[0]} -f {input[1]} -o {output} &> {log}
        """

rule EventMatch:
    input:
        BUILD_DIR + "/102_EventMatch",
        DUMP_DIR + "/101_EventRecon/Run{run_id}_recon.root"
    output:
        DUMP_DIR + "/102_EventMatch/Run{run_id}_matched.root"
    log:
        LOG_DIR + "/102_EventMatch/Run{run_id}_match.log"
    shell:
        """
        mkdir -p {LOG_DIR}
        mkdir -p {DUMP_DIR}/102_EventMatch
        {input[0]} -f {input[1]} -o {output} &> {log}
        """

# rule rootify_all:
#     input:
#         expand(DUMP_DIR + "/102_EventMatch/Run{run_id}_matched.root", run_id=run_ids)
#     output:
#         DUMP_DIR + "/all_matched.done"
#     shell:
#         """
#         touch {output}
#         """