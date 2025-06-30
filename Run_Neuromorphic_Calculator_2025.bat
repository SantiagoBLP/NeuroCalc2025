@echo off
title Neuromorphic Calculator 2025 by Dr. Ing. Santiago Barrionuevo

:: ---------------------------------------------------------------------
:: Ensure that the working directory is this batch file's folder.
:: This lets generate_qe_eels_inputs.py find ..\inputs\scf.in correctly.
:: ---------------------------------------------------------------------
setlocal enableextensions
cd /d %~dp0

set ERRFILE=outputs\errors.log

:: Define QE executables path (assuming "qe" subfolder in the same directory as this script)
set QE=qe

:: ======================================================
:: Set Fullscreen Mode via VBScript
:: ======================================================
set TempVBSFile=%temp%\~tmpSendKeysTemp.vbs
echo Set WshShell = WScript.CreateObject("WScript.Shell") > "%TempVBSFile%"
echo Wscript.Sleep 500 >> "%TempVBSFile%"
echo WshShell.SendKeys "{F11}" >> "%TempVBSFile%"
cscript //nologo "%TempVBSFile%"
del "%TempVBSFile%"

:: -------------------------------
:: 1. Display Welcome Screen
:: -------------------------------
python python_dependencies\welcome.py
cls

:: -------------------------------
:: Input Method Selection Menu
:: -------------------------------
:input_menu
cls
echo ====================================================
echo    Neuromorphic Calculator 2025 - Input Options
echo ====================================================
echo [1] I have an ID for a Materials Project / Let's simulate
echo [2] I don't know my Materials Project ID / I need assistance
echo [3] I'm trying to simulate a new unpublished material and need help
echo ====================================================
set /p option="Enter your choice (1, 2, or 3): "
if "%option%"=="1" goto option1
if "%option%"=="2" goto option2
if "%option%"=="3" goto option3
echo Invalid choice. Try again.
timeout /t 2 > nul
goto input_menu

:option1
cls
echo =============================================================
echo Standard Workflow: Simulation using a Materials Project ID
echo =============================================================
echo.
echo Please enter the Materials Project ID of the crystal you want to simulate:
set /p MATERIAL_ID="Material ID: "
echo Generating Quantum ESPRESSO input files for material %MATERIAL_ID%...
python python_dependencies\generate_qe_inputs.py %MATERIAL_ID%
if %errorlevel% neq 0 (
    echo ERROR: Failed to generate QE inputs. Possibly a parameter was not properly defined.
    echo ERROR: Failed to generate QE inputs. >> "%ERRFILE%"
    pause
    goto calc_fail
)
goto simulation_workflow

:option2
cls
echo =============================================================
echo Assistance: Identifying a Known Material
echo =============================================================
echo Deploying AI infrastructure...one sec!
python python_dependencies\assist_known_material.py
if %errorlevel% neq 0 (
    echo ERROR: Assistance process failed. Possibly a parameter was not properly defined.
    echo ERROR: Assistance process failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)

:: Extract material_id from selected_materials.json
for /f "tokens=2 delims=:, " %%a in ('findstr /i "material_id" selected_materials.json') do set "MATERIAL_ID=%%a"
:: Remove any quotes from the extracted id
set "MATERIAL_ID=%MATERIAL_ID:"=%"

echo Generating QE input files for material %MATERIAL_ID%...
python python_dependencies\generate_qe_inputs.py %MATERIAL_ID%
if %errorlevel% neq 0 (
    echo ERROR: Failed to generate QE inputs. Possibly a parameter was not properly defined.
    echo ERROR: Failed to generate QE inputs. >> "%ERRFILE%"
    pause
    goto calc_fail
)
goto simulation_workflow

:option3
cls
echo =============================================================
echo Expert Workflow: Unpublished Material Simulation
echo =============================================================
echo Deploying Expert AI infrastructure... This will take a moment one sec!
python python_dependencies\assist_unpublished_material.py
if %errorlevel% neq 0 (
    echo ERROR: Expert assistance process failed. Possibly a parameter was not properly defined.
    echo ERROR: Expert assistance process failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)
echo Generating QE input files from ad-hoc parameters...
python python_dependencies\generate_qe_inputs_ad_hoc.py
if %errorlevel% neq 0 (
    echo ERROR: Failed to generate QE inputs. Possibly a parameter was not properly defined.
    echo ERROR: Failed to generate QE inputs. >> "%ERRFILE%"
    pause
    goto calc_fail
)
goto simulation_workflow

:simulation_workflow
cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [*] SCF Calculation (Running...)
echo [-] NSCF Calculation
echo [-] Bands Calculation
echo [-] Bands Post-Processing
echo [-] Band Plotting
echo [-] DOS Calculation
echo [-] DOS Plotting
echo [-] EELS Low Loss Spectra
echo [-] EELS Plotting
echo =============================================================
echo Your simulation is running in the background. You can check the progress in outputs\scf.out.
echo Simulating material %MATERIAL_ID%...

:: Define MPI process count
set NP=6

:: --- SCF Calculation ---
mpiexec -np %NP% %QE%\pw.exe < inputs\scf.in > outputs\scf.out 2>&1
if %errorlevel% neq 0 (
    echo ERROR: SCF calculation failed. Possibly a parameter was not properly defined.
    echo ERROR: SCF calculation failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)

cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [*] NSCF Calculation (Running...)
echo [-] Bands Calculation
echo [-] Bands Post-Processing
echo [-] Band Plotting
echo [-] DOS Calculation
echo [-] DOS Plotting
echo [-] EELS Low Loss Spectra
echo [-] EELS Plotting
echo =============================================================
echo Your NSCF simulation is running in the background. You can check the progress in outputs\nscf.out.
echo Simulating material %MATERIAL_ID%...

:: --- NSCF Calculation ---
mpiexec -np %NP% %QE%\pw.exe < inputs\nscf.in > outputs\nscf.out 2>&1
if %errorlevel% neq 0 (
    echo ERROR: NSCF calculation failed. Possibly a parameter was not properly defined.
    echo ERROR: NSCF calculation failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)

cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [OK] NSCF Calculation (Completed)
echo [*] Bands Calculation (Running...)
echo [-] Bands Post-Processing
echo [-] Band Plotting
echo [-] DOS Calculation
echo [-] DOS Plotting
echo [-] EELS Low Loss Spectra
echo [-] EELS Plotting
echo =============================================================
echo Your BANDS CALCULATION is running in the background. You can check the progress in outputs\bands.out.
echo Bands for material %MATERIAL_ID%...

:: --- Bands Calculation ---
mpiexec -np %NP% %QE%\pw.exe < inputs\bands.in > outputs\bands.out 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Bands calculation failed. Possibly a parameter was not properly defined.
    echo ERROR: Bands calculation failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)

cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [OK] NSCF Calculation (Completed)
echo [OK] Bands Calculation (Completed)
echo [*] Bands Post-Processing (Running...)
echo [-] Band Plotting
echo [-] DOS Calculation
echo [-] DOS Plotting
echo [-] EELS Low Loss Spectra
echo [-] EELS Plotting
echo =============================================================
echo Your BANDS POST-PROCESSING is running in the background. You can check the progress in outputs\bands_pp.out.
echo Bands for material %MATERIAL_ID%...

:: --- Bands Post-Processing ---
mpiexec -np %NP% %QE%\bands.exe < inputs\bands_pp.in > outputs\bands_pp.out 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Band post-processing failed. Possibly a parameter was not properly defined.
    echo ERROR: Band post-processing failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)

cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [OK] NSCF Calculation (Completed)
echo [OK] Bands Calculation (Completed)
echo [OK] Bands Post-Processing (Completed)
echo [*] Band Plotting (Running...)
echo [-] DOS Calculation
echo [-] DOS Plotting
echo [-] EELS Low Loss Spectra
echo [-] EELS Plotting
echo =============================================================
echo Plotting Bands! Material %MATERIAL_ID%...
echo =============================================================

:: --- Band Plotting ---
python python_dependencies\plot_bands.py 2>> "%ERRFILE%"
if %errorlevel% neq 0 (
    echo ERROR: Band plotting failed. Possibly a parameter was not properly defined.
    echo ERROR: Band plotting failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)

cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [OK] NSCF Calculation (Completed)
echo [OK] Bands Calculation (Completed)
echo [OK] Bands Post-Processing (Completed)
echo [OK] Band Plotting (Completed)
echo [*] DOS Calculation (Running...)
echo [-] DOS Plotting
echo [-] EELS Low Loss Spectra
echo [-] EELS Plotting
echo =============================================================
echo Your DOS CALCULATION is running in the background. You can check the progress in outputs\dos.out.
echo Simulating DOS for material %MATERIAL_ID%...

:: --- DOS Calculation ---
mpiexec -np %NP% %QE%\dos.exe < inputs\dos.in > outputs\dos.out 2>&1
if %errorlevel% neq 0 (
    echo ERROR: DOS calculation failed. Possibly a parameter was not properly defined.
    echo ERROR: DOS calculation failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)

cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [OK] NSCF Calculation (Completed)
echo [OK] Bands Calculation (Completed)
echo [OK] Bands Post-Processing (Completed)
echo [OK] Band Plotting (Completed)
echo [OK] DOS Calculation (Completed)
echo [*] DOS Plotting (Running...)
echo [-] EELS Low Loss Spectra
echo [-] EELS Plotting
echo =============================================================
echo Plotting DOS! Material %MATERIAL_ID%...
echo =============================================================

:: --- DOS Plotting ---
python python_dependencies\plot_dos.py 2>> "%ERRFILE%"
if %errorlevel% neq 0 (
    echo ERROR: DOS plotting failed. Possibly a parameter was not properly defined.
    echo ERROR: DOS plotting failed. >> "%ERRFILE%"
    pause
    goto calc_fail
)

cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [OK] NSCF Calculation (Completed)
echo [OK] Bands Calculation (Completed)
echo [OK] Bands Post-Processing (Completed)
echo [OK] Band Plotting (Completed)
echo [OK] DOS Calculation (Completed)
echo [OK] DOS Plotting (Completed)
echo [*] EELS Low Loss Spectra (Running...)
echo [-] EELS Plotting
echo =============================================================
echo Calculating EELS Low Loss Spectra! Material %MATERIAL_ID%...
echo =============================================================

:: --- EELS Calculation and Post-Processing ---
echo Your EELS CALCULATION is running in the background. You can check the progress in outputs\qe_eels_low.out
echo Simulating EELS for material %MATERIAL_ID%...

python python_dependencies\generate_qe_eels_inputs.py
if %errorlevel% neq 0 (
    echo ERROR: Failed to generate QE EELS inputs. Possibly a parameter was not properly defined.
    echo ERROR: Failed to generate QE EELS inputs. >> "%ERRFILE%"
    pause
    goto calc_fail
)

echo Running low-q EELS (turbo_eels.x) on %NP% cores...
mpiexec -np %NP% turbo_eels.exe < inputs\qe_eels_low.in > outputs\qe_eels_low.out 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Low Loss EELS calculation failed!
    echo ERROR: Low Loss EELS calculation failed! >> "%ERRFILE%"
    pause
    goto calc_fail
)

echo Running low-q EELS post-processing (turbo_spectrum.x) on %NP% cores...
mpiexec -np %NP% turbo_spectrum.exe < inputs\qe_eels_low_spectrum.in > outputs\qe_eels_low_spectrum.out 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Low Loss EELS post-processing failed!
    echo ERROR: Low Loss EELS post-processing failed! >> "%ERRFILE%"
    pause
    goto calc_fail
)

cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [OK] NSCF Calculation (Completed)
echo [OK] Bands Calculation (Completed)
echo [OK] Bands Post-Processing (Completed)
echo [OK] Band Plotting (Completed)
echo [OK] DOS Calculation (Completed)
echo [OK] DOS Plotting (Completed)
echo [OK] EELS Low Loss Spectra (Completed)
echo [*] EELS Plotting (Running...)
echo =============================================================
echo Plotting EELS Low Loss Spectra! Material %MATERIAL_ID%...
echo =============================================================
echo Now Plotting EELS Spectra...
python python_dependencies\plot_all_eels.py 2>> "%ERRFILE%"

if %errorlevel% neq 0 (
    echo ERROR: EELS plotting failed!
    echo ERROR: EELS plotting failed! >> "%ERRFILE%"
    pause
    goto calc_fail
)

:: ------------------- Final Simulation Summary -------------------
cls
echo =============================================================
echo                    Simulation Workflow
echo =============================================================
echo [OK] SCF Calculation (Completed)
echo [OK] NSCF Calculation (Completed)
echo [OK] Bands Calculation (Completed)
echo [OK] Bands Post-Processing (Completed)
echo [OK] Band Plotting (Completed)
echo [OK] DOS Calculation (Completed)
echo [OK] DOS Plotting (Completed)
echo [OK] EELS Low-Loss Spectra Calculation (Completed)
echo =============================================================
echo CALCULATION COMPLETE CONGRATS!
echo Data processing finished successfully.
echo =============================================================
echo You can check your results in outputs\.

:: ======================================================
:: 7. AI Analysis of QE Outputs via Perplexity API
:: ------------------------------------------------------
echo Initiating automated scientific analysis of QE outputs...
python python_dependencies\AI_perp.py

:: ======================================================
:: 8. Ask the User to Retry or Exit
:: ======================================================
:retry_or_exit
cls
echo What would you like to do next?
echo [1] Retry with a new simulation workflow
echo [2] Exit the program
set /p choice="Enter your choice (1 or 2): "
if "%choice%"=="1" (
    goto input_menu
) else if "%choice%"=="2" (
    cls
    echo Exiting the program. Keep exploring the quantum world!
    timeout /t 2 > nul
    exit
) else (
    echo Invalid choice! Please enter 1 or 2.
    timeout /t 2 > nul
    goto retry_or_exit
)

:calc_fail
cls
echo ================================================================
echo ERROR: The calculation has failed. Possibly a parameter was not properly defined.
echo Do you want to start over?
echo [1] Yes, start over.
echo [2] No, exit.
echo ================================================================
set /p fail_choice="Enter your choice (1 or 2): "
if "%fail_choice%"=="1" (
   goto input_menu
) else if "%fail_choice%"=="2" (
   echo Exiting the program. Keep exploring the quantum world!
   timeout /t 2 > nul
   exit /b
) else (
   echo Invalid choice!
   timeout /t 2 > nul
   goto calc_fail
)
