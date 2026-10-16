import os, urllib.request, urllib.error

op = urllib.request.build_opener()
op.addheaders = [('User-Agent', 'Mozilla/5.0 (X11; Linux x86_64)')]
urllib.request.install_opener(op)

MW, LOGP, PURCH = 'CDEFGHIJ', 'ABCDE', 'AB'
os.makedirs('data/zinc', exist_ok=True)
got = miss = 0
codes = {}
for m in MW:
    for l in LOGP:
        for p in PURCH:
            name = f"{m}{l}A{p}"
            dst = f"data/zinc/{name}.smi"
            if os.path.exists(dst) and os.path.getsize(dst) > 0:
                got += 1; continue
            try:
                urllib.request.urlretrieve(
                    f"https://files.docking.org/2D/{m}{l}/{name}.smi", dst)
                got += 1
                print(f"  {name}  {os.path.getsize(dst)//1024} KB", flush=True)
            except urllib.error.HTTPError as e:
                miss += 1; codes[e.code] = codes.get(e.code, 0) + 1
            except Exception as e:
                miss += 1; codes[type(e).__name__] = codes.get(type(e).__name__, 0) + 1

print(f"\ngot {got} files, {miss} absent   errors: {codes}")

seen, n = set(), 0
with open('data/zinc_instock.smi', 'w') as out:
    for f in sorted(os.listdir('data/zinc')):
        for line in open(f'data/zinc/{f}'):
            q = line.split()
            if len(q) < 2 or q[0].lower() == 'smiles' or q[0] in seen:
                continue
            seen.add(q[0]); out.write(f"{q[0]} {q[1]}\n"); n += 1
print(f"merged {n} unique molecules -> data/zinc_instock.smi")
