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

echo Installing / verifying dependencies...
python -m pip install -q -r requirements.txt

echo.
echo Choose execution mode:
echo [1] Launch Interactive Streamlit Web App (Recommended)
echo [2] Run Batch CLI Pipeline in Terminal
echo.
set /p choice="Enter choice (1 or 2, default is 1): "

if "%choice%"=="2" (
    echo Running GelPhase CLI Pipeline...
    python run_analysis.py
) else (
    echo Launching Streamlit Web App at http://localhost:8501 ...
    python -m streamlit run app.py
)

echo.
pause
