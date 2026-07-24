# Comprehensive audit of manuscript
import re, os

os.chdir("D:/Researching/NHANES")

with open("manuscript.tex", encoding="utf-8") as f:
    text = f.read()
    lines = text.split("\n")

print("FILE INFO")
print(f"  Lines: {len(lines)}")
print(f"  Chars: {len(text)}")

# Sections
secs = [l.strip() for l in lines if l.strip().startswith("\\"+"section{") and not l.strip().startswith("\\"+"section*")]
print(f"\nSECTIONS: {len(secs)}")
for s in secs:
    print(f"  {s[:90]}")

# Tables
tcounter = 0
for i, l in enumerate(lines):
    if "\\begin{table}" in l:
        tcounter += 1
        for j in range(i, min(i+8, len(lines))):
            if "\\caption{" in lines[j]:
                print(f"  Table {tcounter}: {lines[j][:100]}")
                break
print(f"\nTOTAL TABLES: {tcounter}")

# Citations - use manual parsing
cites = []
for l in lines:
    idx = l.find("\\cite{")
    while idx >= 0:
        end = l.find("}", idx)
        if end >= 0:
            content = l[idx+6:end]
            cites.extend([k.strip() for k in content.split(",")])
        idx = l.find("\\cite{", idx+1)
unique = sorted(set([c for c in cites if c]))
print(f"\nUNIQUE CITATIONS: {len(unique)}")
print(f"  Keys: {', '.join(unique)}")

# Check for TESTCUTOFF
if "TESTCUTOFF" in text:
    print("\n!!! WARNING: TESTCUTOFF marker still in file!")

# Check for placeholders
for p in ["TODO", "FIXME", "XXX", "TBD"]:
    if p in text:
        print(f"WARNING: '{p}' in text!")

# Supplementary listing
print(f"\n=== SUPPLEMENTARY TABLE LISTING ===")
with open("supplementary.tex", encoding="utf-8") as f:
    stext = f.read()
slines = stext.split("\n")
for i, l in enumerate(slines):
    if "\\section*{" in l or "Table S" in l:
        print(f"  {l.strip()[:90]}")

# References.bib cross-check
print(f"\n=== REFERENCES CROSS-CHECK ===")
with open("references.bib", encoding="utf-8") as f:
    bib = f.read()
bib_entries = []
for l in bib.split("\n"):
    if l.startswith("@article{"):
        key = l.split("{")[1].rstrip(",")
        bib_entries.append(key)
bib_set = set(bib_entries)

cited_set = set(unique)
missing = cited_set - bib_set
unused = bib_set - cited_set
if missing:
    print(f"  CITED but NOT in bib: {missing}")
if unused:
    print(f"  In bib but NOT cited: {unused}")
if not missing and not unused:
    print(f"  OK: all {len(unique)} cited keys match bib entries.")

# Compilation
print(f"\n=== COMPILATION LOG ===")
with open("manuscript.log", encoding="utf-8", errors="replace") as f:
    log = f.read()
errs = log.count("! ")
undef = log.count("undefined")
print(f"  Errors (!): {errs}")
print(f"  Undefined refs: {undef}")

PDF = "manuscript.pdf"
if os.path.exists(PDF):
    print(f"  PDF: {os.path.getsize(PDF)/1024:.0f} KB")
PDF2 = "supplementary.pdf"
if os.path.exists(PDF2):
    print(f"  Supp PDF: {os.path.getsize(PDF2)/1024:.0f} KB")

print(f"\n{'='*50}")
print("AUDIT COMPLETE")
