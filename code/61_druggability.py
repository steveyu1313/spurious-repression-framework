# -*- coding: utf-8 -*-
"""Annotate the 20 hub genes and the 4 panel genes against DGIdb
(drug-gene interactions) to establish which of them are actionable."""
import io, json, os, time, urllib.request

base = "08_QC_Revision"
tab = os.path.join(base, "tables")

hubs = [l.strip() for l in io.open(os.path.join(tab, "T17_hub_genes.csv"), encoding="utf-8").read().split("\n")[1:] if l.strip()]
panel = ["AMACR", "PCA3", "HOXC6", "HPN"]
genes = sorted(set(hubs) | set(panel))
print("genes:", len(genes))

url = "https://dgidb.org/api/graphql"
query = """
{ genes(names: %s) {
    nodes { name
      interactions {
        drug { name conceptId }
        interactionScore
        interactionTypes { type }
        sources { sourceDbName }
      } } } }
""" % json.dumps(genes)

body = json.dumps({"query": query}).encode("utf-8")
req = urllib.request.Request(url, data=body, headers={"Content-Type": "application/json"})
raw = urllib.request.urlopen(req, timeout=120).read()
d = json.loads(raw)
nodes = d["data"]["genes"]["nodes"]
print("returned genes:", len(nodes))

rows = []
for nd in nodes:
    g = nd["name"]
    ints = nd.get("interactions") or []
    drugs, srcs, types = {}, set(), set()
    for it in ints:
        dn = (it.get("drug") or {}).get("name")
        if not dn:
            continue
        drugs[dn] = max(drugs.get(dn, 0), it.get("interactionScore") or 0)
        for s in it.get("sources") or []:
            if s.get("sourceDbName"):
                srcs.add(s["sourceDbName"])
        for t in it.get("interactionTypes") or []:
            if t.get("type"):
                types.add(t["type"])
    top = sorted(drugs.items(), key=lambda kv: -kv[1])[:5]
    rows.append(dict(
        gene=g,
        n_interactions=len(ints),
        n_drugs=len(drugs),
        n_sources=len(srcs),
        top_drugs="; ".join("%s (%.2f)" % (k, v) for k, v in top),
        sources="; ".join(sorted(srcs)),
        interaction_types="; ".join(sorted(types)),
        in_hub=int(g in hubs),
        in_panel=int(g in panel),
        hub_rank=(hubs.index(g) + 1) if g in hubs else "",
    ))

rows.sort(key=lambda r: (-r["n_drugs"], -r["n_interactions"]))
out = os.path.join(tab, "T36_druggability_DGIdb.csv")
cols = list(rows[0].keys())
with io.open(out, "w", encoding="utf-8", newline="") as f:
    f.write(",".join(cols) + "\n")
    for r in rows:
        vals = []
        for c in cols:
            v = str(r[c])
            if any(ch in v for ch in [",", '"', ";"]):
                v = '"' + v.replace('"', '""') + '"'
            vals.append(v)
        f.write(",".join(vals) + "\n")

print("\n%-10s %5s %5s %5s  %s" % ("gene", "ints", "drugs", "srcs", "top drug"))
for r in rows:
    print("%-10s %5d %5d %5d  %s" % (r["gene"], r["n_interactions"], r["n_drugs"],
                                     r["n_sources"], r["top_drugs"][:60]))
print("\nwritten:", out)

no_drug = [r["gene"] for r in rows if r["n_drugs"] == 0]
print("no annotated drug:", no_drug)
print("median drugs per gene (annotated):",
      sorted([r["n_drugs"] for r in rows if r["n_drugs"] > 0])[len([r for r in rows if r["n_drugs"] > 0]) // 2])
