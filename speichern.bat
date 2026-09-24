@echo off
rem Doppelklick: Selbsttest -> bei OK committen und auf GitHub pushen.
cd /d "%~dp0"

py selbsttest.py > "%TEMP%\nlfem_selbsttest.txt"
type "%TEMP%\nlfem_selbsttest.txt"

findstr /L /C:"[FEHLER" "%TEMP%\nlfem_selbsttest.txt" >nul
if not errorlevel 1 (
    echo.
    echo Selbsttest hat FEHLER - es wurde NICHTS gespeichert.
    pause
    exit /b 1
)

echo.
git status --short
echo.
set "MSG="
set /p "MSG=Was hast du geaendert? (kurz, dann Enter): "
if "%MSG%"=="" set "MSG=Update"

git add .
git commit -m "%MSG%"
git push

echo.
pause
