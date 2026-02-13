import glob
import json

pede_config_files_eeemcal = glob.glob("data/DESY_2025/config/Pedestal_*.json")
pede_config_files_focal = glob.glob("data/SPS_2024/config/Pedestal_*.json")

analysis_adc_config_files_eeemcal = glob.glob("data/DESY_2025/config/AnalysisADC_*.json")
analysis_adc_config_files_focal = glob.glob("data/SPS_2024/config/AnalysisADC_*.json")

# Extract all samples from filenames
pede_samples_eeemal = [file.split('/')[-1].replace('.json','') for file in pede_config_files_eeemcal]
pede_samples_focal  = [file.split('/')[-1].replace('.json','') for file in pede_config_files_focal]

analysis_adc_samples_eeemal = [file.split('/')[-1].replace('.json','') for file in analysis_adc_config_files_eeemcal]
analysis_adc_samples_focal  = [file.split('/')[-1].replace('.json','') for file in analysis_adc_config_files_focal]

qc_all_files_eeemcal = glob.glob("data/DESY_2025/data/beam/Run*.h2g")
qc_all_files_focal = glob.glob("data/SPS_2024/data/beam/Run*.h2g")

adc_analysis_all_files = glob.glob("dump/102_EventMatch/Run*_matched.root")

qc_all_files_eeemcal = [file.split('/')[-1].replace('.h2g','') for file in qc_all_files_eeemcal]
qc_all_files_focal = [file.split('/')[-1].replace('.h2g','') for file in qc_all_files_focal]
adc_analysis_runs = [file.split('/')[-1].replace('_matched.root','') for file in adc_analysis_all_files]

rule Pedestal_EEEMCal:
    input:
        config="data/DESY_2025/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Run Files"],
        script="build/202_PedestalSub"
    output:
        rootfile="dump/202_PedestalSub/EEEMCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/202_PedestalSub/EEEMCal/{sample}.log"
    shell:
        "build/202_PedestalSub -f {input.config} -o {output.rootfile} -e {params.event} &> {log}"

rule Pedestal_FoCal:
    input:
        config="data/SPS_2024/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Run Files"],
        script="build/202_PedestalSub"
    output:
        rootfile="dump/202_PedestalSub/FoCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/202_PedestalSub/FoCal/{sample}.log"
    shell:
        "build/202_PedestalSub -f {input.config} -o {output.rootfile} -e {params.event} --focal &> {log}"

rule Pedestal_EEEMCal_all:
    input:
        expand("dump/202_PedestalSub/EEEMCal/{sample}.root", sample=pede_samples_eeemal)

rule Pedestal_FoCal_all:
    input:
        expand("dump/202_PedestalSub/FoCal/{sample}.root", sample=pede_samples_focal)

rule QC_EEEMCal:
    input:
        script="build/301_QC",
        rootfile="dump/102_EventMatch/EEEMCal/{sample}_matched.root"
    output:
        rootfile="dump/301_QC/EEEMCal/{sample}_qc.root"
    log:
        "logs/301_QC/EEEMCal/{sample}.log"
    shell:
        "build/301_QC -f {input.rootfile} -o {output.rootfile} &> {log}"

rule QC_FoCal:
    input:
        script="build/301_QC",
        rootfile="dump/102_EventMatch/FoCal/{sample}_matched.root"
    output:
        rootfile="dump/301_QC/FoCal/{sample}_qc.root"
    log:
        "logs/301_QC/FoCal/{sample}.log"
    shell:
        "build/301_QC -f {input.rootfile} -o {output.rootfile} &> {log}"

rule QC_EEEMCal_all:
    input:
        expand("dump/301_QC/EEEMCal/{sample}_qc.root", sample=qc_all_files_eeemcal)

rule QC_FoCal_all:
    input:
        expand("dump/301_QC/FoCal/{sample}_qc.root", sample=qc_all_files_focal)

rule AnalysisADC_EEEMCal:
    input:
        config="data/DESY_2025/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Run Files"],
        script="build/203_AnalysisADC",
        pede=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Pedestal Subtraction"],
    output:
        rootfile="dump/203_AnalysisADC/EEEMCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/203_AnalysisADC/EEEMCal/{sample}.log"
    shell:
        "build/203_AnalysisADC -f {input.config} -o {output.rootfile} -e {params.event} -p {input.pede} &> {log}"

rule AnalysisADC_FoCal:
    input:
        config="data/SPS_2024/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Run Files"],
        script="build/203_AnalysisADC",
        pede=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Pedestal Subtraction"],
    output:
        rootfile="dump/203_AnalysisADC/FoCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/203_AnalysisADC/FoCal/{sample}.log"
    shell:
        "build/203_AnalysisADC -f {input.config} -o {output.rootfile} -e {params.event} -p {input.pede} --focal &> {log}"

rule AnalysisADC_EEEMCal_all:
    input:
        expand("dump/203_AnalysisADC/EEEMCal/{sample}.root", sample=analysis_adc_samples_eeemal)

rule AnalysisADC_FoCal_all:
    input:
        expand("dump/203_AnalysisADC/FoCal/{sample}.root", sample=analysis_adc_samples_focal)

rule ADC_Analysis:
    input:
        script="build/303_ADC_Analysis",
        rootfile="dump/102_EventMatch/{run}_matched.root"
    output:
        rootfile="dump/303_ADC_Analysis/{run}.root"
    log:
        "logs/303_ADC_Analysis/{run}.log"
    shell:
        "{input.script} -f {input.rootfile} -o {output.rootfile} &> {log}"

rule ADC_Analysis_all:
    input:
        expand("dump/303_ADC_Analysis/{run}.root", run=adc_analysis_runs)

rule TemplateFit_FoCal:
    input:
        script="build/205_TemplateFit",
        inputfile="dump/102_EventMatch/FoCal/{sample}_matched.root",
        csvfile="data/SPS_2024/config/PhaseScan_2V5_20250207_152421_Chn77.csv",
        timewalkfile="dump/204_TimeWalk/FoCal/{sample}_timewalk.json"
    output:
        rootfile="dump/205_TemplateFit/FoCal/{sample}_fit.root",
        jsonfile="dump/205_TemplateFit/FoCal/{sample}_fit.json"
    # number of CPU cores to use
    threads: 16
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/205_TemplateFit/FoCal/{sample}.log"
    shell:
        "{input.script} -f {input.inputfile} -o {output.rootfile} -c {input.csvfile} -e {params.event} --focal -t {input.timewalkfile} &> {log}"

rule TemplateFit_EEEMCal:
    input:
        script="build/205_TemplateFit",
        inputfile="dump/102_EventMatch/EEEMCal/{sample}_matched.root",
        csvfile="data/DESY_2025/config/PhaseScan_2V5_20250207_152421_Chn77.csv",
        timewalkfile="dump/204_TimeWalk/EEEMCal/{sample}_timewalk.json"
    output:
        rootfile="dump/205_TemplateFit/EEEMCal/{sample}_fit.root",
        jsonfile="dump/205_TemplateFit/EEEMCal/{sample}_fit.json"
    # number of CPU cores to use
    threads: 16
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/205_TemplateFit/EEEMCal/{sample}.log"
    shell:
        "{input.script} -f {input.inputfile} -o {output.rootfile} -c {input.csvfile} -e {params.event} -t {input.timewalkfile} &> {log}"

rule EventFit_FoCal:
    input:
        script="build/208_EventFit",
        inputfile="dump/102_EventMatch/FoCal/{sample}_matched.root",
        csvfile="data/SPS_2024/config/PhaseScan_2V5_20250207_152421_Chn77.csv",
        timewalkfile="dump/204_TimeWalk/FoCal/Run770_timewalk.json"
    output:
        rootfile="dump/208_EventFit/FoCal/{sample}_fit.root"
    # number of CPU cores to use
    threads: 30
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/208_EventFit/FoCal/{sample}.log"
    shell:
        "{input.script} -f {input.inputfile} -o {output.rootfile} -c {input.csvfile} -e {params.event} --focal -t {input.timewalkfile} &> {log}"

rule Timewalk_FoCal:
    input:
        script="build/204_TimeWalk",
        inputfile="dump/102_EventMatch/FoCal/{sample}_matched.root",
        csvfile="data/SPS_2024/config/PhaseScan_2V5_20250207_152421_Chn77.csv"
    output:
        rootfile="dump/204_TimeWalk/FoCal/{sample}_timewalk.root",
        jsonfile="dump/204_TimeWalk/FoCal/{sample}_timewalk.json"
    params:
        event_timewalk = lambda wildcards: config.get("event_timewalk", -1)
    log:
        "logs/204_TimeWalk/FoCal/{sample}_timewalk.log"
    shell:
        "{input.script} -f {input.inputfile} -o {output.rootfile} -c {input.csvfile} -e {params.event_timewalk} --focal --timewalk &> {log}"

rule Timewalk_EEEMCal:
    input:
        script="build/204_TimeWalk",
        inputfile="dump/102_EventMatch/EEEMCal/{sample}_matched.root",
        csvfile="data/DESY_2025/config/PhaseScan_2V5_20250207_152421_Chn77.csv"
    output:
        rootfile="dump/204_TimeWalk/EEEMCal/{sample}_timewalk.root",
        jsonfile="dump/204_TimeWalk/EEEMCal/{sample}_timewalk.json"
    params:
        event_timewalk = lambda wildcards: config.get("event_timewalk", -1)
    log:
        "logs/204_TimeWalk/EEEMCal/{sample}_timewalk.log"
    shell:
        "{input.script} -f {input.inputfile} -o {output.rootfile} -c {input.csvfile} -e {params.event_timewalk} --timewalk &> {log}"

rule AnalysisFit_FoCal:
    input:
        config="data/SPS_2024/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Run Files"],
        script="build/206_AnalysisFit",
        pede=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Pedestal Subtraction"],
        timewalk=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["TimeWalk JSON"],
        fitting=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Fitting JSON"],
        csvfile="data/SPS_2024/config/PhaseScan_2V5_20250207_152421_Chn77.csv"
    output:
        rootfile="dump/206_AnalysisFit/FoCal/{sample}.root"
    params:
        event_analysis = lambda wildcards: config.get("event_analysis", -1)
    log:
        "logs/206_AnalysisFit/FoCal/{sample}.log"
    shell:
        "build/206_AnalysisFit -f {input.config} -o {output.rootfile} -e {params.event_analysis} -p {input.pede} -c {input.csvfile} --focal &> {log}"

rule AnalysisFit_EEEMCal:
    input:
        config="data/DESY_2025/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Run Files"],
        script="build/206_AnalysisFit",
        pede=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Pedestal Subtraction"],
        timewalk=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["TimeWalk JSON"],
        fitting=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Fitting JSON"],
        csvfile="data/DESY_2025/config/PhaseScan_2V5_20250207_152421_Chn77.csv"
    output:
        rootfile="dump/206_AnalysisFit/EEEMCal/{sample}.root"
    params:
        event_analysis = lambda wildcards: config.get("event_analysis", -1)
    log:
        "logs/206_AnalysisFit/EEEMCal/{sample}.log"
    shell:
        "build/206_AnalysisFit -f {input.config} -o {output.rootfile} -e {params.event_analysis} -p {input.pede} -c {input.csvfile} &> {log}"