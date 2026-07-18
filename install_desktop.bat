@echo off
setlocal
cd /d "%~dp0"

set "PYTHON="
set "PYTHON_ARGS="
if exist "%USERPROFILE%\miniconda3\python.exe" set "PYTHON=%USERPROFILE%\miniconda3\python.exe"
if not defined PYTHON if exist "%USERPROFILE%\anaconda3\python.exe" set "PYTHON=%USERPROFILE%\anaconda3\python.exe"
if not defined PYTHON (
  where python >nul 2>nul
  if not errorlevel 1 set "PYTHON=python"
)
if not defined PYTHON (
  where py >nul 2>nul
  if not errorlevel 1 set "PYTHON=py"
  if not errorlevel 1 set "PYTHON_ARGS=-3"
)

if not defined PYTHON (
  echo Python was not found. Install Python 3.11 or newer, then run this installer again.
  pause
  exit /b 1
)

"%PYTHON%" %PYTHON_ARGS% -m pip install -r requirements-desktop.txt
if errorlevel 1 (
  echo Installation failed. Check Python and network access.
  pause
  exit /b 1
)
echo SpatialTX Studio Desktop dependencies are ready.
pause
