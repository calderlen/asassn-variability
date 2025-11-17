# Changelog

## 17-11-2025
- Added a CLI entry point directly to `calder/df_filter_naive.py`, eliminating the need for the separate `script_filter_naive.py`.
- Removed `calder/script_filter_naive.py` to avoid duplicated logic.
- Added build_manifest to search for asas-sn ID's in lc_cal subdirs BEFORE doing reproduction check, so reproduce_candidates only needs to check against a csv file instead of searching through all mag_bins each time.