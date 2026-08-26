# Shihai's Laser Injection Data Analysis

## Introduction

This data analysis is based on the KCU105-based data analysis for the laser injection experiment conducted in the FoCal Lab at CERN.

#### Test setup

The test setup is based on the H2GCROC3D ASICs, Xilinx KCU105 FPGA board, and a laser injection system.
The laser signal is **NOT** synchronized with the clock of the ASICs, and the laser itself is difused by a piece of tape, which distributes the light over a larger area of the sensor but **NOT** evenly.

<img src="doc/LaserSetup_Photo.png" alt="Test setup" height="320" />
<img src="doc/LaserSetup_Scheme.png" alt="Test setup schematic" height="320" />

## Data Processing

The subfolders of the project is organized as follows:

```
.
├── doc/                # Documentation and images
├── data/               # Raw data files
├── src/                # Source code for libraries
├── include/            # Header files for libraries
├── logs/               # Log files from data processing
├── config/             # Configuration files
├── build/              # Build files for the project
├── scripts/            # Scripts for data processing and analysis
├── dump/               # Generated files from data processing
├── workflow/           # Snakemake workflow rule files for automation
└── readme.md           # This file
````

### ROOT conversion

The raw data files are converted to ROOT format using the `101_EventRecon.cxx`, `102_EventMatch.cxx`, and `103_Rootifier.cxx` scripts. These are the same scripts used in the previous KCU-based beam tests.
