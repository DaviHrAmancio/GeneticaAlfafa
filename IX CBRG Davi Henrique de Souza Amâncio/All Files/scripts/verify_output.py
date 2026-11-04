#!/usr/bin/env python3
"""Confere os resultados-chave do snapshot contra os CSVs e contra os números do resumo submetido."""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
snap = json.loads((ROOT / "config_snapshot.json").read_text())
k = snap["metricas_chave"]
chk = pd.read_csv(ROOT / "outputs" / "csv" / "verificacao_csv_vs_bruto.csv").iloc[0]

checks = [
    ("65 observações", k["n_observacoes"] == 65),
    ("CSV idêntico ao LAID bruto", chk.n_pareadas == 65 and chk.max_dif_abs_LAI == 0),
    ("P12 AGZG1219 = 0,600 (resumo submetido)", abs(k["P12_AGZG1219_csv"] - 0.600) < 1e-3),
    ("P55 AGZG1501 = 0,333 (resumo submetido)", abs(k["P55_AGZG1501_csv"] - 0.333) < 1e-3),
    ("Frobenius = 1,0650 (resumo submetido)", abs(k["Frobenius_treino_csv"] - 1.0650) < 1e-4),
]
ok = True
for name, passed in checks:
    print(("OK   " if passed else "FALHA"), name)
    ok &= bool(passed)

for f in ["observacoes_lai", "transicoes", "matrizes_transicao_long", "distancias_frobenius", "ciclos_rebrota", "contexto_alfafa_fontes"]:
    p = ROOT / "outputs" / "csv" / f"{f}.csv"
    exists = p.exists() and p.stat().st_size > 0
    print(("OK   " if exists else "FALHA"), f"{f}.csv")
    ok &= exists
sys.exit(0 if ok else 1)
