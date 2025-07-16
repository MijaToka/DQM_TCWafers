# TCWafers creation
This repository handles the creation of TC wafer templates for the [HGCal DQM](https://gitlab.cern.ch/hgcal-dpg/hgcal-comm).

The TC are made via a jupyter notebook `generateTCwafer.ipynb` where the TC's are defined as rombuses in the sequence given by the `WaferCellMapTraces.txt` file. The ordering of these cells was chosen by visual inspection thanks to the `generateTCvisualization.ipynb` which holds the scripts to read this input file and process it to display the module so that the value of the cell corresponds to the value of the "Global TC count", ie the index asigned to the TC by the Cell indexer.
This file also contains python sketches of the [Cell](https://github.com/MijaToka/cmssw/blob/dev/trggmapping_hgcal_up2date/CondFormats/HGCalObjects/interface/HGCalMappingCellIndexerTrigger.h) and [Module](https://github.com/MijaToka/cmssw/blob/dev/trggmapping_hgcal_up2date/CondFormats/HGCalObjects/interface/HGCalMappingModuleIndexerTrigger.h) indexer for the TC and an [ESProducer](https://github.com/MijaToka/cmssw/blob/dev/trggmapping_hgcal_up2date/Geometry/HGCalMapping/plugins/HGCalMappingTriggerESProducer.cc) sketch that later was included into CMSSW.

## Instalation
```bash
git clone https://github.com/MijaToka/DQM_TCWafers.git
```
## Python dependencies
- ROOT
- NumPy
- Pandas
