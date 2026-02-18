# Universal-Plasmid-Maker


This repository implements a simple **Universal Plasmid Maker** tool for the coding assignment.  
Given:

- An input DNA sequence in FASTA format (`Input.Fa`, e.g. a genome or a plasmid like `pUC19.fa`)
- A design file (`Design.txt`) describing the desired multiple cloning site (MCS) and antibiotic markers
- A marker dictionary file (`markers.tab`) listing antibiotic markers and their DNA sequences

the tool outputs a new plasmid sequence (`Output.Fa`) that is compatible with the input organism.

## Features

- Reads a single‑record DNA FASTA file.
- Heuristically identifies an **origin of replication (ORI)** as the most AT‑rich window.
- Automatically appends a **default replication gene cassette** to satisfy the “genes necessary for replication” requirement.
- Parses a **design file** with:
  - `Multiple_Cloning_SiteX, RestrictionEnzymeName`
  - `Antibiotic_markerY, AntibioticName`
- Builds an **MCS** from common restriction sites (EcoRI, BamHI, HindIII, etc.) using an internal dictionary.
- Loads **antibiotic marker sequences** from `markers.tab`, skipping missing markers with a warning instead of failing.
- Can **remove specified restriction sites** (e.g. EcoRI) from the final plasmid sequence by mutating the recognition motif.
- Includes a **test** that validates EcoRI removal for the `pUC19` example.

## File layout

- `plasmid_maker.py` – main library and command‑line tool.
- `test_pUC19.py` – test that uses `pUC19.fa` and `Design_pUC19.txt`.
- `pUC19.fa` – example input plasmid sequence (provided by instructor).
- `Design_pUC19.txt` – example design file for pUC19 (provided by instructor).
- `markers.tab` – marker dictionary file (provided by instructor).

