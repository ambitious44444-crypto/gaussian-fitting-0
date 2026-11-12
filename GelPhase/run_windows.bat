@echo off
chcp 65001 > nul
echo =====================================================================
echo  🧬 GelPhase: Automated Gel Lane Profile Deconvolution
echo =====================================================================
echo Checking Python environment...

python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python is not installed or not in PATH!
    echo Please install Python 3.9+ from https://www.python.org/
    pause
    exit /b
)

echo Installing dependencies if needed...
python -m pip install -q -r requirements.txt

echo Running GelPhase Pipeline...
python run_analysis.py

echo.
echo =====================================================================
echo  Analysis Complete! Check the 'results/' folder for CSVs and plots.
echo =====================================================================
pause
