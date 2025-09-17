import subprocess
import sys

required_packages = [
    "numpy",
    "matplotlib",
    "scipy",
    "pymatgen",
    "mp-api",
    "colorama",
    "openai",
]

# Optional for Windows
if sys.platform.startswith("win"):
    required_packages.append("windows-curses")

def install_package(package):
    try:
        print(f"🔍 Checking {package}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f"✅ {package} installed.")
    except subprocess.CalledProcessError:
        print(f"❌ Failed to install {package}")

def main():
    for pkg in required_packages:
        install_package(pkg)

if __name__ == "__main__":
    main()
