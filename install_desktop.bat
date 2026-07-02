@echo off
setlocal
cd /d "%~dp0"

set "PYTHON=python"
if exist "%USERPROFILE%\anaconda3\python.exe" set "PYTHON=%USERPROFILE%\anaconda3\python.exe"
if exist "%USERPROFILE%\miniconda3\python.exe" set "PYTHON=%USERPROFILE%\miniconda3\python.exe"

"%PYTHON%" -m pip install -r requirements-desktop.txt
if errorlevel 1 (
  echo Installation failed. Check Python and network access.
  pause
  exit /b 1
)
echo SpatialTX Studio Desktop dependencies are ready.
pause
