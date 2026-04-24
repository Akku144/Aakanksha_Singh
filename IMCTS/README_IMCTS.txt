IMCTS integration notes

- This folder was added to the NEW copy:
  /Users/aakankshasingh/Downloads/Preliminary-Mission-Design-Tool-Public-main-IMCTS

- Original folder was not modified:
  /Users/aakankshasingh/Downloads/Preliminary-Mission-Design-Tool-Public-main

- Run from MATLAB:
  cd('/Users/aakankshasingh/Downloads/Preliminary-Mission-Design-Tool-Public-main-IMCTS')
  out = run_combined_pmdt_imcts(40, 10000, 2025, true, 'fuel');
  out = run_combined_pmdt_imcts(40, 10000, 2025, true, 'time');
  out = run_combined_pmdt_imcts(40, 10000, 2025, true, 'weighted', 0.02);

- IMCTS code and data loader are in:
  /Users/aakankshasingh/Downloads/Preliminary-Mission-Design-Tool-Public-main-IMCTS/IMCTS

- Notes:
  - `run_combined_pmdt_imcts` uses PMDT debris catalog (`Data/DebrisData.mat`).
  - PMDT `getPosition.m` now supports all debris entries in DebrisData, not only 3 IDs.
