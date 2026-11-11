# cms-raw-csv

Writes the raw readout of CMS events as plain CSV files: one directory per event, one file per subsystem. Every row
is one channel that was read out, with its address, its raw values (every time sample) and its position in cm. The
subsystems are pixel, strip, ECAL, preshower, HCAL, DT, CSC, RPC and GEM, plus the size of every FED's data block.

The input is the RAW data tier. The standard CMSSW `RawToDigi` sequence unpacks it, and the EDAnalyzer
`EventCsvDumper` (`plugins/EventCsvDumper.cc`) writes the files.

[`csv_map.html`](csv_map.html) shows in 3D which file covers which part of the detector. Open it in a browser; it loads
three.js from a CDN.

## Run it from scratch

All you need is `/cvmfs/cms.cern.ch`, network access to the CMS conditions database (see
[Conditions](#conditions)) and about 2 GB of memory. Tested with CMSSW_14_0_25 on 2024 data.

```bash
# 0. CMSSW_14_0_25 is built for EL9 (RHEL, Alma or Rocky Linux 9). On any other OS, open an EL9 container first and
#    run everything below inside it. It needs apptainer or singularity on the host, and it sees your home and the
#    top-level directory you start it from; add -B /some/dir for anything else, e.g. where RAW.root is.
/cvmfs/cms.cern.ch/common/cmssw-el9

# 1. A CMSSW release area
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_14_0_25
cd CMSSW_14_0_25/src
cmsenv

# 2. This package, built (a minute or two)
git clone https://github.com/sgnoohc/cms-raw-csv.git RawViz/EventCsv
scram b -j 4

# 3. Dump the first event of RAW.root, from the directory RAW.root is in
cd /path/to/directory/with/RAW.root
cmsRun $CMSSW_BASE/src/RawViz/EventCsv/test/dump_csv_cfg.py inputFiles=file:RAW.root
ls event_csv/
```

This writes `event_csv/run<run>_evt<event>/`. The first run on a machine takes a few minutes while the conditions
are fetched; after that a run takes about 30 seconds, most of it spent loading conditions and geometry.

Later, in a new shell (inside `cmssw-el9` again if you needed it), set up and go straight to step 3:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /path/to/CMSSW_14_0_25/src && cmsenv
```

## Choosing events

```bash
edmFileUtil -e RAW.root | head -30                                  # list the run, lumi section and event numbers
cmsRun $CMSSW_BASE/src/RawViz/EventCsv/test/dump_csv_cfg.py inputFiles=file:RAW.root maxEvents=5              # the first 5
cmsRun $CMSSW_BASE/src/RawViz/EventCsv/test/dump_csv_cfg.py inputFiles=file:RAW.root events=383631:994851190  # chosen ones
```

| Option | Default | Meaning |
|---|---|---|
| `inputFiles` | `file:RAW.root` | RAW files, comma-separated. `file:` marks a local path; `root://...` works too. |
| `maxEvents` | `1` | Dump the first this many events (`-1` = all). Ignored when `events` is given. |
| `events` | | Dump these events, `run:event`, comma-separated. |
| `outputDir` | `event_csv` | Each event goes into `<outputDir>/run<run>_evt<event>/`. |
| `era` | `Run3_2024` | Data-taking era. Use `Run3` for 2022 and `Run3_2023` for 2023 data. |
| `globalTag` | `auto:run3_data_prompt` | Conditions, which cover every Run 3 run so far. |
| `threads` | `4` | Threads for the unpacking. |

Only 2024 data has been tested. For 2025 data or Run 2 data, use a CMSSW release from that period.

## Conditions

The unpackers need cabling maps, and the dumper needs the detector geometry. `cmsRun` reads both from the CMS
Frontier conditions database, through the squid proxies listed in
`/cvmfs/cms.cern.ch/SITECONF/local/JobConfig/site-local-config.xml`. The machine's cvmfs configuration
(`CMS_LOCAL_SITE`) decides where `local` points. At CMS sites and on lxplus it already works.

If `cmsRun` stops with an error about `site-local-config.xml` or Frontier, tell it which site's configuration to use.
The squids of that site must accept connections from your machine:

```bash
ls /cvmfs/cms.cern.ch/SITECONF/                             # the sites
export SITECONFIG_PATH=/cvmfs/cms.cern.ch/SITECONF/T2_CH_CERN   # e.g. on a CERN machine
```

## What you get

```
event_csv/run383631_evt994851190/
  pixel.csv strip.csv ecal_barrel.csv ecal_endcap.csv ecal_preshower.csv hcal_hbhe.csv hcal_hf.csv hcal_ho.csv
  dt.csv csc_strips.csv csc_wires.csv csc_comparators.csv rpc.csv gem.csv fed.csv
  extent/      the outline of every channel, row for row with the CSV of the same name
  inventory/   every readout element of the detector, read out or not
```

- Each CSV opens with `#` comment lines: what the file holds, then one line per column. The row of column names
  comes after them. Read a file with `pandas.read_csv(path, comment='#')`.
- Every row is one channel that was read out in this event. Channels that were not read out are absent; most were
  dropped by zero suppression.
  - ECAL and the CSC strips are read out in blocks around activity, so many of their rows hold only pedestal.
  - HF is read out in full every event.
- All values are raw codes as stored, before calibration.
- Every file except `fed.csv` ends in `x_cm, y_cm, z_cm`: the channel's position in global CMS coordinates, in cm.
  x points toward the centre of the LHC, y up, and z along the beam.
  - Pixel, strip, CSC, RPC and GEM: the centre of the pixel, strip or wire group.
  - DT: the midpoint of the wire.
  - ECAL and HCAL: the cell's reference point from the CMSSW geometry, on the cell's front face.
- `extent/<file>` gives, for every row of `<file>`, the four corners of the channel's sensitive area, in order
  around its edge (`x0_cm` ... `z3_cm`).
  - A strip is outlined over its whole length.
  - A DT hit is outlined as its 4.2 cm drift cell along the wire.
  - A CSC wire group is outlined as its band of wires.
  - A calorimeter cell is outlined by its face toward the collision point.
- `inventory/` lists every readout element with its address and its size in channels:
  - the pixel and strip modules, with their centres and axes;
  - DT and CSC layers, RPC rolls and GEM eta partitions;
  - the ECAL endcap, preshower and HCAL cells. The ECAL barrel is a full 170 × 360 grid and is not listed.
- HCAL calibration and monitoring channels are left out.

| File | Columns | Values |
|---|---|---|
| `pixel.csv` | part (BPix/FPix), layer (BPix layer 1–4, FPix disk 1–3), detid, row, col, adc | adc: 8-bit, 0–255. One number per pixel; no time samples. |
| `strip.csv` | part (TIB/TID/TOB/TEC), layer, stereo (0/1), detid, strip, adc | adc: 0–253 is linear; 254 means 254–1022; 255 means saturated. One number per strip. |
| `ecal_barrel.csv` | ieta, iphi, adc0–adc9, gain0–gain9 | 10 samples, 25 ns apart, 12-bit ADC (0–4095). The gain is written as the factor 12, 6 or 1. The pedestal is about 200 counts at gain 12; an in-time pulse peaks near sample 5. |
| `ecal_endcap.csv` | zside, ix, iy, adc0–adc9, gain0–gain9 | Same as the barrel. |
| `ecal_preshower.csv` | zside, plane, six, siy, strip, adc0–adc2 | 3 samples. Pedestal-subtracted, so values can be negative. |
| `hcal_hbhe.csv` | part (HB/HE), ieta, iphi, depth (1–7), soi, adc0–adc7, tdc0–tdc7 | 8 samples, 25 ns apart. soi is the index of the sample from the triggered crossing. adc is an 8-bit nonlinear QIE11 code; turning it into charge needs a lookup table from the conditions database. HE tdc: 6-bit in 0.5 ns steps, 0–49 valid, 62/63 mean no valid time. HB tdc: a 2-bit code, where 3 means none. |
| `hcal_hf.csv` | ieta, iphi, depth (1–4), soi, adc0–adc2, tdc_rise0–2, tdc_fall0–2 | 3 samples. adc is an 8-bit QIE10 code. tdc_rise: 6-bit, 0.5 ns steps, 62/63 mean no valid time. tdc_fall: 5-bit. |
| `hcal_ho.csv` | ieta, iphi, depth (4), soi, adc0–adc9 | 10 samples. 7-bit QIE8 code, 0–127. |
| `dt.csv` | wheel, station, sector, superlayer, layer, wire, hit, tdc_counts, time_ns | No amplitude. time_ns = tdc_counts × 25/32 ns, still including the trigger latency. A wire can have several hits, indexed by `hit`. |
| `csc_strips.csv` | zside (±1), station, ring, chamber, layer, strip, adc0–adc7 | 8 samples, 50 ns apart, 12-bit. The pedestal is about 500–1000 and differs per strip. ME1/1 is split: ring 1 = ME1/1b and ring 4 = ME1/1a. |
| `csc_wires.csv` | zside, station, ring, chamber, layer, wiregroup, timebins_on | timebins_on lists the 25 ns bins (0–15) that fired, space-separated. |
| `csc_comparators.csv` | zside, station, ring, chamber, layer, strip, half (0/1 = left/right half-strip), timebins_on | Binary hit plus the time bins that fired; used by the trigger. ME1/1 comparators are all in ring 1: strips 65–112 are ME1/1a strips 1–48, and their positions are on ME1/1a. |
| `rpc.csv` | region (0 = barrel, ±1 = endcap), ring (wheel in the barrel), station, sector, layer, subsector, roll, strip, bx | Hit or not, plus the bunch crossing relative to the triggered one. |
| `gem.csv` | region (±1), station, ring, chamber, layer, roll (eta partition), strip, bx | Hit or not, plus the bunch crossing. |
| `fed.csv` | fed_id, subsystem, bytes | The size of each FED's raw data block, header and trailer included. FEDs missing from `FEDNumbering.h` are labelled `other`. |

The comment lines at the top of each file repeat these descriptions in more detail.
