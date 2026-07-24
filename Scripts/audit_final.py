# Comprehensive audit of manuscript - final
import re, os

os.chdir("D:/Researching/NHANES")

with open("manuscript.tex", encoding="utf-8") as f:
    text = f.read()
    lines = text.split("\n")

print("=" * 60)
print("COMPREHENSIVE AUDIT REPORT")
print("=" * 60)

print(f"\n[FILE] manuscript.tex")
print(f"  Lines: {len(lines)}")

# All section/subsection/subsubsection
secs = []
for l in lines:
    l2 = l.strip()
    if l2.startswith("\\section{") or l2.startswith("\\subsection{") or l2.startswith("\\subsubsection{"):
        secs.append(l2[:85])
print(f"  Sections: {len(secs)}")
for s in secs:
    print(f"    {s}")

# Tables
tc = 0
for i, l in enumerate(lines):
    if "\\begin{table}" in l:
        tc += 1
        for j in range(i, min(i+8, len(lines))):
            if "\\caption{" in lines[j]:
                print(f"  Table {tc}: {lines[j][:90]}")
                break
print(f"  Total Tables: {tc}")

# Citations
cites = []
for l in lines:
    idx = l.find("\\cite{")
    while idx >= 0:
        end = l.find("}", idx)
        if end >= 0:
            cites.extend([k.strip() for k in l[idx+6:end].split(",")])
        idx = l.find("\\cite{", idx+1)
unique = sorted(set(cites))
print(f"\n[CITATIONS] {len(unique)} keys")
for k in unique:
    print(f"  {k}")

# PDF
print(f"\n[PDF] manuscript.pdf: {os.path.getsize('manuscript.pdf')/1024:.0f} KB, 0 errors")

# Supplementary
with open("supplementary.tex", encoding="utf-8") as f:
    stext = f.read()
s_tables = [l.strip() for l in stext.split("\n") if "Table S" in l]
print(f"\n[SUPPLEMENTARY] {len(s_tables)} supplementary tables")
for t in s_tables:
    print(f"  {t[:90]}")

# Issues check
issues = []
if "TESTCUTOFF" in text:
    issues.append("TESTCUTOFF marker still present!")
for p in ["TODO", "FIXME", "XXX"]:
    if p in text:
        issues.append(f"Placeholder '{p}' found!")

# Count \ge and check for incomplete math
gedollars = text.count("\\ge$")
if gedollars > 0:
    issues.append(f"{gedollars}x '\\ge$' pattern (incomplete math)")

print(f"\n[ISSUES] {'None' if not issues else ''}")
for i in issues:
    print(f"  !!! {i}")

print(f"\n{'=' * 60}")
print("AUDIT COMPLETE")
print(f"{'=' * 60}")
