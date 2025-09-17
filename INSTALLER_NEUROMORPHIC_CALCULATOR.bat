:: === Step 0: Check and install Python dependencies ===
echo ====================================================
echo Checking Python dependencies...
echo ====================================================
python python_dependencies\check_and_install_requirements.py
if %errorlevel% neq 0 (
    echo ❌ Failed to verify/install dependencies. Check pip configuration.
    pause
    exit /b
)
