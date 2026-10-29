# fetch_hpa_images.py — Human Protein Atlas IHC for panel & core genes
import re, json, time, urllib.request, pathlib

OUT = pathlib.Path(r'D:/yuproject/前列腺癌/08_QC_Revision/data/external/hpa')
OUT.mkdir(parents=True, exist_ok=True)
HDRS = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}

GENES = {
    'AMACR': 'ENSG00000242110',
    'HOXC6': 'ENSG00000197757',
    'HPN':   'ENSG00000105707',
    'MYC':   'ENSG00000136997',
    'EZH2':  'ENSG00000106462',
    'FASN':  'ENSG00000169710',
}

def get(url, binary=False, retries=2):
    for a in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=HDRS)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read() if binary else r.read().decode('utf-8', 'ignore')
        except Exception as e:
            if a == retries: raise
            time.sleep(2)

manifest = {}
for gene, eid in GENES.items():
    manifest[gene] = {}
    # --- prostate adenocarcinoma IHC ---
    try:
        h = get(f'https://www.proteinatlas.org/{eid}-{gene}/cancer/prostate+cancer')
        urls = re.findall(r'images\.proteinatlas\.org/(\d+)/([0-9A-Za-z_]+)\.jpg', h)
        can = [f'{a}/{b}' for a, b in urls if not b.endswith('_thumb')]
        can = sorted(set(can))
        manifest[gene]['cancer'] = can[:4]
        for i, p in enumerate(can[:4]):
            data = get(f'https://images.proteinatlas.org/{p}.jpg', binary=True)
            (OUT / f'{gene}_cancer_{i}.jpg').write_bytes(data)
            print(f'{gene} cancer {i}: {p} ({len(data)//1024} KB)')
    except Exception as e:
        print(f'{gene} CANCER FAIL: {e}')
        manifest[gene]['cancer'] = []
    # --- normal prostate IHC ---
    try:
        h = get(f'https://www.proteinatlas.org/{eid}-{gene}/tissue/prostate')
        urls = re.findall(r'images\.proteinatlas\.org/(\d+)/([0-9A-Za-z_]+)_selected_medium\.jpg', h)
        nor = sorted(set(f'{a}/{b}' for a, b in urls if 'rna' not in b))
        manifest[gene]['normal'] = nor[:4]
        for i, p in enumerate(nor[:4]):
            data = get(f'https://images.proteinatlas.org/{p}_selected_medium.jpg', binary=True)
            (OUT / f'{gene}_normal_{i}.jpg').write_bytes(data)
            print(f'{gene} normal {i}: {p} ({len(data)//1024} KB)')
    except Exception as e:
        print(f'{gene} NORMAL FAIL: {e}')
        manifest[gene]['normal'] = []
    time.sleep(1)

# PCA3 check (lncRNA, expect no IHC)
try:
    s = json.loads(get('https://www.proteinatlas.org/search/PCA3?format=json'))
    manifest['PCA3_HPA'] = [r.get('Gene') for r in s][:3]
except Exception as e:
    manifest['PCA3_HPA'] = f'query failed: {e}'
print('PCA3 in HPA:', manifest['PCA3_HPA'])

(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
print('manifest written')
