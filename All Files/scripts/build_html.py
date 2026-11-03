#!/usr/bin/env python3
"""Gera html/painel_alfafa.html (arquivo único) a partir de outputs/csv/*.csv e dos .MOW brutos."""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pipeline_alfafa_lai as P  # noqa: E402

CSV = ROOT / "outputs" / "csv"
obs = pd.read_csv(CSV / "observacoes_lai.csv")
mt = pd.read_csv(CSV / "matrizes_transicao_long.csv")
perm = pd.read_csv(CSV / "teste_permutacao_entre_experimentos.csv")
stab = pd.read_csv(CSV / "estabilidade_temporal_2017.csv")
sens = pd.read_csv(CSV / "sensibilidade_regra_de_estados.csv")
cyc = pd.read_csv(CSV / "resumo_ciclos.csv")

mows = {}
for g in P.GENOS:
    t = P.read_dssat_table(P.RAW / f"{g}.MOW")
    mows[g] = [P.dssat_date(d).isoformat() for d in t["DATE"].dropna() if P.dssat_date(d).year in (2016, 2017)]

data = {
    "obs": obs[["experimento", "data", "LAI", "estado_csv", "estado_resumo", "amostra"]].to_dict("records"),
    "mows": mows,
    "mat": mt[(mt.subconjunto == "todas")][["regra_estados", "experimento", "amostra", "de", "para", "contagem", "n_linha", "prob", "ic95_inf", "ic95_sup"]]
        .round(5).to_dict("records"),
    "perm": perm.round(4).to_dict("records"),
    "stab": stab.round(4).to_dict("records"),
    "sens": sens.round(4).fillna("").to_dict("records"),
    "cyc": cyc.round(4).fillna("").to_dict("records"),
    "estados": P.STATE_NAMES,
    "edges": P.EDGES,
}
html = (ROOT / "scripts" / "painel_template.html").read_text(encoding="utf-8")
out = ROOT / "html" / "painel_alfafa.html"
out.parent.mkdir(exist_ok=True)
out.write_text(html.replace("/*__DATA__*/null", json.dumps(data, ensure_ascii=False, default=str)), encoding="utf-8")
print("OK", out, f"{out.stat().st_size/1024:.0f} KB")
