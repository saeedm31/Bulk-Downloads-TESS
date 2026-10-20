# 🛰️ Bulk-Downloads-TESS

**Download TESS light curves for a filtered list of stars in one go, instead of clicking through the archive.**

Built during my PhD to pull light curves for every bright star observed in TESS sectors 1–13.

## How it works

1. **Combine target lists:** reads the official TESS target lists for sectors 1–13 (`store/all_targets_S0XX_v1.txt`) and removes stars that appear in several sectors.
2. **Filter:** keeps the stars you want, by default those brighter than **Tmag 9**. You can also filter by RA and Dec.
3. **Build a download script:** matches those TIC IDs against MAST's [bulk-download scripts](https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html) and writes a single `final_download.sh` with only the light curves you need.
4. **Download:**

```bash
python create_download_query.py
chmod +x Data_storage/final_download.sh
./Data_storage/final_download.sh
```

`main_tess_search.py` goes further. For each selected star it searches light curves with `lightkurve`, queries stellar parameters from MAST and VizieR, and computes variability measures (weighted RMS, bolometric corrections, trigonometric log g).

**Requirements:** `pandas`, `numpy`, `scipy`, `astropy`, `astroquery`, `lightkurve`.

> Written for sectors 1–13. For later sectors, add their target lists and `tesscurl_sector_N_lc.sh` files to `store/` and extend the sector range.

---

Part of my astronomy background: [Saeed Hojjatpanah](https://github.com/saeedm31), PhD in Astronomy (exoplanets, TESS).
