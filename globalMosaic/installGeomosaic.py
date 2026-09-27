"""
Console-script wrapper around install/installGeomosaic.

The real work is a shell script, because it is doing git and make: keeping it as a shell
script means it can also be run on a bare machine that has no Python environment yet, which
is the usual case on a fresh instance. This module only lets `pip install -e .` put it on
PATH as `installGeomosaic`, and forwards every argument unchanged.
"""
import os
import sys
from pathlib import Path


def scriptPath():
    """Locate install/installGeomosaic, which sits beside the package, not inside it."""
    return Path(__file__).resolve().parent.parent / "install" / "installGeomosaic"


def main():
    script = scriptPath()
    if not script.is_file():
        # Only reachable from a wheel/non-editable install, where the sibling install/
        # directory is not shipped. Say so rather than failing obscurely.
        sys.stderr.write(
            f"installGeomosaic: cannot find {script}\n"
            "  This wrapper resolves the script relative to the source tree, so it needs an\n"
            "  editable install (pip install -e .). From a source checkout, run\n"
            "  install/installGeomosaic directly instead.\n")
        return 1
    if not os.access(script, os.X_OK):
        os.chmod(script, 0o755)
    # execv so the shell script owns the tty and its exit status is ours: a caller using
    # subprocess.run(check=True) must see a build or capability failure.
    os.execv("/bin/bash", ["bash", str(script)] + sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main())
