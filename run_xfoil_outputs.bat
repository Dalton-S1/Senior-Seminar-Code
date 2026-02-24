@echo off
setlocal EnableExtensions EnableDelayedExpansion

echo =====================================
echo Running XFOIL scripts 1-650 (60s timeout each)
echo =====================================

for /L %%i in (1,1,650) do (
  if exist "output_%%i.txt" (
    echo ---------------------------------
    echo Running output_%%i.txt

    powershell -NoProfile -ExecutionPolicy Bypass -Command ^
      "$p = Start-Process -FilePath 'xfoil.exe' -NoNewWindow -PassThru -RedirectStandardInput 'output_%%i.txt';" ^
      "if (Wait-Process -Id $p.Id -Timeout 60 -ErrorAction SilentlyContinue) {" ^
      "  Write-Host 'XFOIL finished normally.'" ^
      "} else {" ^
      "  Write-Host 'XFOIL timed out at 60s — terminating.';" ^
      "  Stop-Process -Id $p.Id -Force" ^
      "}"

    timeout /t 1 >nul
  )
)

echo =====================================
echo All XFOIL runs completed.
echo =====================================
pause
