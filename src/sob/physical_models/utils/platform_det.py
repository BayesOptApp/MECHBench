r""" 
THIS IS A MODULE TO AID DETECTING THE PLATFORM OF THE SYSTEM AND SETTING UP 
THE PATH TO OPENRADIOSS, AS WELL AS OTHER RUNNER OPTIONS.
"""

# Module Properties
__author__ = "Ivan Olarte Rodriguez"


import platform
import zipfile
import os
import requests
from pathlib import Path
import stat

ALLOWED_PLATFORMS = {'Windows', 'Linux'}

# Cache platform once
_SYSTEM = platform.system()


def platform_detection() -> str:
    if _SYSTEM == 'Darwin':
        return 'MacOS'
    return _SYSTEM


def is_allowed_platform() -> bool:
    return _SYSTEM in ALLOWED_PLATFORMS


def raise_if_not_allowed_platform():
    if not is_allowed_platform():
        raise ValueError(
            f"Unsupported platform: {_SYSTEM}. Allowed: {ALLOWED_PLATFORMS}"
        )


def _get_download_info():
    if _SYSTEM == 'Windows':
        return (
            "https://github.com/OpenRadioss/OpenRadioss/releases/download/latest-20260319/OpenRadioss_win64.zip",
            "win"
        )
    elif _SYSTEM == 'Linux':
        return (
            "https://github.com/OpenRadioss/OpenRadioss/releases/download/latest-20260319/OpenRadioss_linux64.zip",
            "linux"
        )
    else:
        raise ValueError(f"Unsupported platform: {_SYSTEM}")


def download_zip_openradioss() -> str:
    raise_if_not_allowed_platform()

    url, system_lower = _get_download_info()

    home_dir = Path.cwd()
    zip_path = home_dir / f"OpenRadioss_{system_lower}64.zip"
    extract_path = home_dir / f"OpenRadioss_{system_lower}64"


    # Stream download (memory efficient)
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(zip_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    # Extract
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_path)
    
    make_executable(extract_path)

    zip_path.unlink(missing_ok=True)

    return str(extract_path / "OpenRadioss")


def make_executable(path: Path):
    for root, _, files in os.walk(path):
        for f in files:
            file_path = Path(root) / f
            file_path.chmod(file_path.stat().st_mode | stat.S_IEXEC)