import subprocess
import os
import sys

VARNA = os.path.join(os.path.dirname(__file__), "../VARNAcmd.jar")

# Waiting for VARNA to be released in a package


def main():
    cmd = ["java", "-jar", VARNA] + sys.argv[1:]
    subprocess.run(cmd)


if __name__ == "__main__":
    main()
