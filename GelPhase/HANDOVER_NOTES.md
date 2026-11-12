# 🤝 AI Handover & Project Continuation Notes (GelPhase)

**Project Name**: GelPhase (Automated Gel Lane Profile Deconvolution & Saturation Threshold Analyzer)  
**Primary Goal**: Open-source computational tool for BIOMOD/iGEM biomolecular design research.  
**Author / User**: 洪有志  
**Last Updated**: 2026-11-13  

---

## 📌 1. Critical User Requirement & Recent Changes
The user specifically requested:
> **"只保留允許使用者匯入的功能（功能 B：拖拉上傳 CSV）。主要原因是希望簡化程式，以確保能夠運作為優先。"**
> *(Only keep the user CSV import feature to simplify the codebase and prioritize ensuring it works reliably on Streamlit Cloud).*

### What Was Just Done:
1. **Completely Simplified `app.py`**:
   - Removed all dependencies on local static folders (`results/`, `data/`, etc.).
   - The app now operates **100% in-memory / on-the-fly** via `st.file_uploader`.
   - Added `sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))` at the top of `app.py` to prevent `ModuleNotFoundError: No module named 'gelphase'` on Linux containers.
2. **`requirements.txt`**:
   - Contains: `streamlit>=1.30.0`, `numpy>=1.24.0`, `scipy>=1.10.0`, `matplotlib>=3.7.0`, `pandas>=2.0.0`.
3. **Workflow of the Simplified App**:
   - User uploads 1 or more ImageJ CSVs (Column 0: pixel $x$, Column 1: intensity $y$).
   - App auto-detects or prompts for concentrations.
   - User clicks **"🚀 Run GelPhase Deconvolution & Threshold Fitting"**.
   - Output: Parameter table (with CSV download), Saturation Threshold plot (with 300 DPI PNG download), and individual multi-Gaussian fits.

---

## 🛠️ 2. Streamlit Cloud Debugging Guide for Incoming AI
If Streamlit Cloud (`share.streamlit.io`) still throws an error:

1. **Check the Streamlit Logs**:
   - Click "Manage App" in bottom-right corner of Streamlit Cloud to see the exact traceback.
2. **Main File Path Setting**:
   - If the GitHub repo root contains `app.py`, Main File Path = `app.py`.
   - If the GitHub repo root is the parent folder, Main File Path = `GelPhase/app.py`.
3. **Git Sync**:
   - Ensure the user has pushed the newly simplified `app.py` to GitHub:
     ```bash
     git add app.py requirements.txt
     git commit -m "Simplify app to upload-only mode"
     git push
     ```
   - In Streamlit Cloud, click "Clear cache and deploy" or "Reboot app".

---

## 📂 3. Directory Layout
```text
GelPhase/
├── app.py                 # Pure interactive web app (upload-only, 0 file dependencies)
├── requirements.txt       # Dependencies including streamlit
├── HANDOVER_NOTES.md      # This file
├── README.md              # Project documentation
├── run_windows.bat        # Local Windows launcher
├── gelphase/              # Core algorithm library
│   ├── __init__.py
│   ├── models.py          # Gaussian & broken-line models
│   ├── deconvolution.py   # Multi-peak deconvolution + AICc
│   ├── threshold.py       # Piecewise linear regression
│   └── visualizer.py      # Matplotlib publication figures
└── results/               # Precomputed benchmark plots (kept locally)
```
