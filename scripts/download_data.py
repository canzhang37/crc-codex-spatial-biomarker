# single-cell table + supplementary files
import hashlib
import subprocess
import sys
import zipfile
from pathlib import Path

CSV = "CRC_clusters_neighborhoods_markers.csv"
CSV_URL = ("https://data.mendeley.com/public-files/datasets/mpjzbtfgfr/files/"
           "c24351b3-76d7-444f-9edf-0246356b0c78/file_downloaded")
CSV_SHA256 = "416cc3926a7a900ce3b22a33be699535ca35ae85fca4585a8ba5dad6d7f3c677"
SUPP_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7479520/supplementaryFiles"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def curl(url, path):
    subprocess.run(["curl", "-s", "-S", "-L", "-C", "-", "--retry", "5", "-o", str(path), url], check=True)


def main(dest):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    p = dest / CSV
    if not p.exists() or sha256(p) != CSV_SHA256:
        curl(CSV_URL, p)
    assert sha256(p) == CSV_SHA256, "sha256 mismatch for the single-cell table"
    print(CSV, "ok")

    z = dest / "PMC7479520_supplementary.zip"
    if not z.exists():
        curl(SUPP_URL, z)
    with zipfile.ZipFile(z) as zf:
        print("supplementary files:", zf.namelist())
        zf.extractall(dest / "supplementary")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/raw")
