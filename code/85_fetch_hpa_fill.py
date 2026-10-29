# fetch_hpa_fill.py — fill missing normal-prostate IHC + retry EZH2 cancer
import re, json, time, urllib.request, pathlib

OUT = pathlib.Path(r'D:/yuproject/前列腺癌/08_QC_Revision/data/external/hpa')
HDRS = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
GENES = {
    'AMACR': 'ENSG00000242110', 'HOXC6': 'ENSG00000197757', 'HPN': 'ENSG00000105707',
    'MYC': 'ENSG00000136997', 'EZH2': 'ENSG00000106462', 'FASN': 'ENSG00000169710',
}

def get(url, binary=False, retries=4):
    last = None
    for a in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=HDRS)
            with urllib.request.urlopen(req, timeout=90) as r:
                data = r.read()
            if binary and len(data) < 20000:
                raise IOError(f'too small: {len(data)}')
            return data if binary else data.decode('utf-8', 'ignore')
        except Exception as e:
            last = e; time.sleep(2)
    raise last

man = json.load(open(OUT / 'manifest.json'))

for gene, eid in GENES.items():
    # normal prostate
    try:
        h = get(f'https://www.proteinatlas.org/{eid}-{gene}/tissue/prostate')
        urls = re.findall(r'images\.proteinatlas\.org/(\d+)/([0-9A-Za-z]+_[0-9A-Za-z]+_[0-9A-Za-z_]+?)\.jpg', h)
        cand = [f'{a}/{b}' for a, b in urls
                if 'rna' not in b and 'selected' not in b and 'thumb' not in b]
        cand = sorted(set(cand))
        man[gene]['normal'] = cand[:4]
        for i, p in enumerate(cand[:4]):
            f = OUT / f'{gene}_normal_{i}.jpg'
            if not f.exists() or f.stat().st_size < 20000:
                data = get(f'https://images.proteinatlas.org/{p}_medium.jpg', binary=True)
                f.write_bytes(data)
                print(f'{gene} normal {i}: {p} ({len(data)//1024} KB)')
    except Exception as e:
        print(f'{gene} NORMAL FAIL: {e}')
    time.sleep(1)

# retry EZH2 cancer (2nd and 3rd images)
try:
    h = get('https://www.proteinatlas.org/ENSG00000106462-EZH2/cancer/prostate+cancer')
    urls = re.findall(r'images\.proteinatlas\.org/(\d+)/([0-9A-Za-z_]+)\.jpg', h)
    can = sorted(set(f'{a}/{b}' for a, b in urls
                     if 'rna' not in b and 'selected' not in b and not b.endswith('_thumb')))
    man['EZH2']['cancer'] = can[:4]
    for i, p in enumerate(can[:4]):
        f = OUT / f'EZH2_cancer_{i}.jpg'
        if not f.exists() or f.stat().st_size < 20000:
            data = get(f'https://images.proteinatlas.org/{p}.jpg', binary=True)
            f.write_bytes(data)
        print(f'EZH2 cancer {i}: {p} ({(OUT / f"EZH2_cancer_{i}.jpg").stat().st_size//1024} KB)')
except Exception as e:
    print(f'EZH2 CANCER FAIL: {e}')

# HOXC6 cancer check
h = get('https://www.proteinatlas.org/ENSG00000197757-HOXC6/cancer/prostate+cancer')
print('HOXC6 page mentions Prostate cancer:', h.count('Prostate cancer'), '| imid count:', len(re.findall(r'imid_\d+', h)))
can = sorted(set(re.findall(r'images\.proteinatlas\.org/(\d+/[0-9A-Za-z_]+)\.jpg', h)))
print('HOXC6 cancer image candidates:', can[:6])

(OUT / 'manifest.json').write_text(json.dumps(man, indent=2))
print('done')
