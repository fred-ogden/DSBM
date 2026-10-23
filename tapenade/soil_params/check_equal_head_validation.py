#!/usr/bin/env python3
"""FLO 10/26: Numerical acceptance checks for DSBM Tapenade experiments."""
import math
import re
import sys
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def numbers(text):
    return [float(x) for x in re.findall(r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?", text)]


def full(path):
    text = Path(path).read_text()
    require("Consistency check passed" in text, "primal/tangent consistency missing")
    require("E10a ADJOINT versus TANGENT" in text, "E10a missing")
    require("E10b ADJOINT versus CENTRAL FINITE DIFFERENCES" in text, "E10b missing")
    a = text.split("E10a ADJOINT versus TANGENT", 1)[1].split("E10b ADJOINT", 1)[0]
    rows = re.findall(r"^dcost/d\S+\s+([+-]?[\d.eE+-]+)\s+([+-]?[\d.eE+-]+)\s+([\d.eE+-]+)", a, re.M)
    require(len(rows) == 8, f"expected 8 E10a parameter comparisons, found {len(rows)}")
    require(all(math.isfinite(float(row[2])) and float(row[2]) < 1e-9 for row in rows), "E10a tangent/adjoint mismatch")
    dots = re.findall(r"^\s+v[1-4]:.*?rel diff\s+([\d.eE+-]+)", a, re.M)
    require(len(dots) == 8, f"expected 8 E10a dot products, found {len(dots)}")
    require(all(math.isfinite(float(x)) and float(x) < 1e-9 for x in dots), "E10a dot-product mismatch")
    b = text.split("E10b ADJOINT versus CENTRAL FINITE DIFFERENCES", 1)[1]
    blocks = re.findall(r"dcost/d\S+\s+adjoint[^\n]*\n(.*?)(?=\ndcost/d|\nLarge steps:|\Z)", b, re.S)
    require(len(blocks) == 4, f"expected four E10b parameter blocks, found {len(blocks)}")
    for i, block in enumerate(blocks):
        matches = re.findall(r"^\s+1e-0[2-6]\s+[+-]?[\d.eE+-]+\s+([\d.eE+-]+)\s*$", block, re.M)
        require(len(matches) == 5, f"E10b parameter {i+1}: expected five FD steps")
        errors = [float(x) for x in matches]
        require(all(math.isfinite(x) for x in errors), "nonfinite FD error")
        require(min(errors) < 1e-5, f"E10b parameter {i+1}: no converged FD step")
    print("PASS: E10a tangent/adjoint and dot products; E10b finite differences")


def active(tangent_path, adjoint_path):
    t = Path(tangent_path).read_text()
    a = Path(adjoint_path).read_text()
    transfer = re.search(r"Limited transfer:\s*([\d.eE+-]+)", t)
    derivative = re.search(r"AD derivative:\s*([\d.eE+-]+)", t)
    require(transfer and derivative, "missing active-limiter tangent output")
    v = float(transfer.group(1))
    d = float(derivative.group(1))
    require(0 < v < 0.015, "limiter did not activate")
    require(math.isfinite(d) and abs(d) > 1e-10, "invalid tangent derivative")
    # Clapp-Hornberger suction and center-to-center elevation difference.
    upper = 0.380 - v / 0.10
    lower = 0.200 + v / 0.30
    psi_up = 0.07275 * (upper / 0.439) ** (-5.5)
    psi_lo = 0.07275 * (lower / 0.439) ** (-5.5)
    require(abs(psi_up - psi_lo - 0.20) < 1e-8, "equal-head residual too large")
    steps = re.findall(r"^\s+1e-0[2-6]\s+[+-]?[\d.eE+-]+\s+([\d.eE+-]+)\s*$", t, re.M)
    require(len(steps) == 5, "missing tangent finite-difference table")
    require(min(float(x) for x in steps) < 1e-7, "tangent/FD mismatch")
    adj = re.search(r"Adjoint derivative:\s*([\d.eE+-]+)", a)
    require(adj, "missing adjoint derivative")
    db = float(adj.group(1))
    require(math.isfinite(db) and abs(db - d) / abs(d) < 1e-10, "tangent/adjoint mismatch")
    require("PASS: active-limiter adjoint derivative" in a, "adjoint driver did not pass")
    print(f"PASS: active limiter; head residual {psi_up-psi_lo-0.20:.3e} m")
    print(f"PASS: tangent/FD and tangent/adjoint; derivative {d:.15e} m")


if __name__ == "__main__":
    try:
        if len(sys.argv) == 3 and sys.argv[1] == "full":
            full(sys.argv[2])
        elif len(sys.argv) == 4 and sys.argv[1] == "active":
            active(sys.argv[2], sys.argv[3])
        else:
            raise ValueError("usage: check_equal_head_validation.py full RESULTS | active TANGENT ADJOINT")
    except (ValueError, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
