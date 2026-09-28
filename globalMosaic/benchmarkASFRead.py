"""Console-script wrapper: the benchmark itself lives in bench/benchmarkASFRead.py."""
import os
import sys
from pathlib import Path


def main():
    script = Path(__file__).resolve().parent.parent / "bench" / "benchmarkASFRead.py"
    if not script.is_file():
        sys.stderr.write("benchmarkASFRead: cannot find %s (needs an editable install)\n" % script)
        return 1
    os.execv(sys.executable, [sys.executable, str(script)] + sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main())
