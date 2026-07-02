@echo off
setlocal
cd /d "%~dp0"

if exist "%USERPROFILE%\miniconda3\python.exe" (
  "%USERPROFILE%\miniconda3\python.exe" desktop_app.py
  if errorlevel 1 goto :launch_error
  goto :eof
)

if exist "%USERPROFILE%\anaconda3\python.exe" (
  "%USERPROFILE%\anaconda3\python.exe" desktop_app.py
  if errorlevel 1 goto :launch_error
  goto :eof
)

where python >nul 2>nul
if %errorlevel%==0 (
  python desktop_app.py
  if errorlevel 1 goto :launch_error
  goto :eof
)

echo Python was not found. Install Python 3.11+ and run install_desktop.bat first.
pause
exit /b 1

:launch_error
echo.
echo SpatialTX Studio Desktop could not start.
echo Run install_desktop.bat, then try again. The error details are shown above.
pause
exit /b 1
