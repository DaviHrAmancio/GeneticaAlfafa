#!/usr/bin/env python3
"""
Pipeline reprodutível — Modelagem markoviana do Índice de Área Foliar (LAI) em alfafa.
IX Congresso Brasileiro de Recursos Genéticos (CBRG).

Fluxo:
  1. Lê os arquivos brutos DSSAT (.ALT = observado, .MOW = cortes, .ALX = experimento).
  2. Confere o CSV "preparado pelo ChatGPT" contra o LAID bruto.
  3. Discretiza o LAI em 5 estados por DUAS regras (limiares do resumo; coluna `estado` do CSV).
  4. Estima matrizes de transição (treino 2016 / teste 2017), incerteza e testes.
  5. Analisa os ciclos de rebrota entre cortes (variáveis fisiológicas, sem discretização).
  6. Gera CSVs, figuras, snapshot de configuração e fontes para o NotebookLM.

Uso:  python pipeline_alfafa_lai.py
Saídas: outputs/csv/*.csv, figures/*.png, config_snapshot.json, notebooklm_sources/
"""
from __future__ import annotations

import hashlib
import json
import platform
import shutil
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "outputs" / "csv"
FIG = ROOT / "figures"
NLM = ROOT / "notebooklm_sources"

SEED = 20261103
N_PERM = 100000
N_BOOT = 100000
GENOS = ["AGZG1219", "AGZG1501"]
EDGES = [1.0, 2.5, 3.8, 5.0]  # limiares do resumo (intervalos fechados à esquerda)
STATE_NAMES = {
    1: "Pós-corte",
    2: "Rebrota inicial",
    3: "Crescimento intermediário",
    4: "Fechamento de dossel",
    5: "Cobertura plena",
}
SCHEMES = {"resumo": "estado_resumo", "csv": "estado_csv"}
SCHEME_LABEL = {
    "resumo": "Limiares do resumo aplicados ao LAI bruto",
    "csv": "Coluna `estado` do CSV (ChatGPT) — reproduz o resumo",
}
SPLITS = ["treino", "teste"]


# ----------------------------------------------------------------------------- leitura DSSAT
def dssat_date(n: int) -> date:
    n = int(n)
    if n >= 1_000_000:  # AAAADDD
        year, doy = n // 1000, n % 1000
    else:  # AADDD
        yy, doy = n // 1000, n % 1000
        year = 2000 + yy if yy < 50 else 1900 + yy
    return date(year, 1, 1) + timedelta(days=doy - 1)


def read_dssat_table(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="latin-1").splitlines()
    h = next(i for i, l in enumerate(lines) if l.startswith("@TRNO"))
    cols = lines[h].replace("@", "").split()
    rows = []
    for l in lines[h + 1:]:
        if not l.strip() or l.startswith(("!", "*")):
            continue
        rows.append(l.split())
    df = pd.DataFrame(rows, columns=cols[: len(rows[0])] if rows else cols)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df.replace(-99, np.nan)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ----------------------------------------------------------------------------- observações
def discretize(lai: float) -> int:
    return int(np.digitize(lai, EDGES, right=False)) + 1


def build_observations() -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    obs_all, mows, cwad_all = [], {}, []
    for g in GENOS:
        alt = read_dssat_table(RAW / f"{g}.ALT")
        mow = read_dssat_table(RAW / f"{g}.MOW")
        mows[g] = sorted(dssat_date(d) for d in mow["DATE"].dropna())
        alt["data"] = alt["DATE"].map(dssat_date)
        # biomassa (CWAD) para produtividade por corte
        cw = alt.loc[alt["CWAD"].notna(), ["data", "CWAD"]].copy()
        cw["experimento"] = g
        cwad_all.append(cw)
        o = alt.loc[alt["LAID"].notna(), ["DATE", "data", "LAID"] + (["CHTD"] if "CHTD" in alt else [])].copy()
        o = o.rename(columns={"LAID": "LAI"})
        o["experimento"] = g
        o = o.sort_values("data").reset_index(drop=True)
        o["year"] = o["data"].map(lambda d: d.year)
        o["amostra"] = np.where(o["year"] == 2016, "treino", "teste")
        # ciclo de corte de cada observação
        cut_ref, pre = [], []
        for i, r in o.iterrows():
            cuts_le = [c for c in mows[g] if c <= r["data"]]
            last = cuts_le[-1] if cuts_le else None
            is_pre = False
            if last is not None and last == r["data"]:
                nxt = o.loc[i + 1, "LAI"] if i + 1 < len(o) else np.nan
                # observação no dia do corte: é pré-corte se a próxima observação cai
                is_pre = bool(nxt < r["LAI"]) if not np.isnan(nxt) else False
                if is_pre:
                    last = cuts_le[-2]
            cut_ref.append(last)
            pre.append(is_pre)
        o["corte_ref"] = cut_ref
        o["obs_pre_corte_no_dia"] = pre
        o["dias_desde_corte"] = [(d - c).days for d, c in zip(o["data"], o["corte_ref"])]
        o["ciclo"] = [f"{g}|{c:%Y-%m-%d}|{y}" for g, c, y in zip(o["experimento"], o["corte_ref"], o["year"])]
        obs_all.append(o)
    obs = pd.concat(obs_all, ignore_index=True)
    obs["estado_resumo"] = obs["LAI"].map(discretize)
    return obs, mows, pd.concat(cwad_all, ignore_index=True)


def check_against_csv(obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    csv = pd.read_csv(RAW / "Alfafa_LAI_Markov_treino_teste.csv", encoding="utf-8-sig")
    csv["data"] = pd.to_datetime(csv["data"]).dt.date
    m = obs.merge(
        csv[["data", "experimento", "LAI", "estado", "amostra"]],
        on=["data", "experimento"], how="outer", suffixes=("", "_csv"), indicator=True,
    )
    both = m[m["_merge"] == "both"]
    chk = pd.DataFrame([{
        "n_obs_bruto": len(obs), "n_obs_csv": len(csv), "n_pareadas": len(both),
        "so_no_bruto": int((m["_merge"] == "left_only").sum()),
        "so_no_csv": int((m["_merge"] == "right_only").sum()),
        "max_dif_abs_LAI": float((both["LAI"] - both["LAI_csv"]).abs().max()),
        "amostra_divergente": int((both["amostra"] != both["amostra_csv"]).sum()),
        "estado_csv_diverge_dos_limiares": int((both["estado"] != both["estado_resumo"]).sum()),
    }])
    obs = obs.merge(csv[["data", "experimento", "estado"]].rename(columns={"estado": "estado_csv"}),
                    on=["data", "experimento"], how="left")
    obs["estado_difere"] = obs["estado_csv"] != obs["estado_resumo"]
    return obs, chk


# ----------------------------------------------------------------------------- transições
def build_transitions(obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (g, s), d in obs.groupby(["experimento", "amostra"]):
        d = d.sort_values("data").reset_index(drop=True)
        for i in range(len(d) - 1):
            a, b = d.loc[i], d.loc[i + 1]
            rows.append({
                "experimento": g, "amostra": s, "data_de": a["data"], "data_para": b["data"],
                "intervalo_dias": (b["data"] - a["data"]).days,
                "LAI_de": a["LAI"], "LAI_para": b["LAI"],
                "estado_resumo_de": a["estado_resumo"], "estado_resumo_para": b["estado_resumo"],
                "estado_csv_de": a["estado_csv"], "estado_csv_para": b["estado_csv"],
                "transicao_com_corte": a["ciclo"] != b["ciclo"], "ciclo_origem": a["ciclo"],
            })
    return pd.DataFrame(rows)


def counts_matrix(t: pd.DataFrame, scheme: str) -> np.ndarray:
    C = np.zeros((5, 5), dtype=int)
    for a, b in zip(t[f"estado_{scheme}_de"], t[f"estado_{scheme}_para"]):
        C[int(a) - 1, int(b) - 1] += 1
    return C


def to_prob(C: np.ndarray) -> np.ndarray:
    n = C.sum(1, keepdims=True)
    return np.divide(C, n, out=np.zeros(C.shape), where=n > 0)


def frob(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.linalg.norm(A - B))


def clopper_pearson(k: int, n: int, alpha=0.05):
    if n == 0:
        return (np.nan, np.nan)
    lo = 0.0 if k == 0 else stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return float(lo), float(hi)


def smooth(C, alpha=0.5):
    S = C + alpha
    return S / S.sum(1, keepdims=True)


def mean_logscore(P: np.ndarray, C: np.ndarray) -> float:
    n = C.sum()
    return float((C * np.log(P)).sum() / n) if n else np.nan


def chain_props(C: np.ndarray):
    n = C.sum(1)
    if (n == 0).any():
        return None
    P = C / n[:, None]
    ncomp, _ = connected_components(P > 0, directed=True, connection="strong")
    irr = ncomp == 1
    pi = np.full(5, np.nan)
    mfpt = np.full(5, np.nan)
    if irr:
        w, v = np.linalg.eig(P.T)
        k = np.argmin(np.abs(w - 1))
        pi = np.real(v[:, k]); pi = pi / pi.sum()
        Q = P[:4, :4]
        m = np.linalg.solve(np.eye(4) - Q, np.ones(4))
        mfpt[:4] = m
        mfpt[4] = 1 + P[4, :4] @ m  # tempo médio de retorno ao estado 5
    return {"P": P, "irredutivel": bool(irr), "pi": pi, "mfpt5": mfpt,
            "permanencia": np.diag(P)}


# ----------------------------------------------------------------------------- análises
def matrices_tables(trans: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    store, rows = {}, []
    for scheme in SCHEMES:
        for g in GENOS:
            for s in SPLITS:
                for sub in ["todas", "sem_corte"]:
                    t = trans[(trans.experimento == g) & (trans.amostra == s)]
                    if sub == "sem_corte":
                        t = t[~t.transicao_com_corte]
                    C = counts_matrix(t, scheme)
                    store[(scheme, g, s, sub)] = C
                    P = to_prob(C)
                    for i in range(5):
                        for j in range(5):
                            lo, hi = clopper_pearson(int(C[i, j]), int(C[i].sum()))
                            rows.append({
                                "regra_estados": scheme, "experimento": g, "amostra": s, "subconjunto": sub,
                                "de": i + 1, "para": j + 1, "contagem": int(C[i, j]),
                                "n_linha": int(C[i].sum()), "prob": round(P[i, j], 6),
                                "ic95_inf": lo, "ic95_sup": hi,
                            })
    return pd.DataFrame(rows), store


def frobenius_tables(store, rng):
    rows, perm_null = [], {}
    # 1) entre genótipos (treino, teste) + permutação por ciclo (treino)
    for scheme in SCHEMES:
        for s in SPLITS:
            F = frob(to_prob(store[(scheme, GENOS[0], s, "todas")]), to_prob(store[(scheme, GENOS[1], s, "todas")]))
            rows.append({"tipo": "entre_experimentos", "regra_estados": scheme, "comparacao": f"{GENOS[0]} vs {GENOS[1]} ({s})",
                         "frobenius": F, "p_permutacao": np.nan, "nota": ""})
        for g in GENOS:
            F = frob(to_prob(store[(scheme, g, "treino", "todas")]), to_prob(store[(scheme, g, "teste", "todas")]))
            rows.append({"tipo": "treino_vs_teste", "regra_estados": scheme, "comparacao": f"{g}: 2016 vs 2017",
                         "frobenius": F, "p_permutacao": np.nan, "nota": "linhas vazias contam como zero"})
    return rows


def permutation_test(trans, scheme, rng):
    """Permuta os rótulos de experimento entre CICLOS de rebrota (treino)."""
    t = trans[trans.amostra == "treino"]
    units = []
    for cyc, d in t.groupby("ciclo_origem"):
        units.append((cyc.split("|")[0], counts_matrix(d, scheme)))
    labels = np.array([u[0] == GENOS[0] for u in units])
    U = np.stack([u[1] for u in units])
    n1 = int(labels.sum())

    def stat(mask):
        A, B = U[mask].sum(0), U[~mask].sum(0)
        return frob(to_prob(A), to_prob(B))

    obs = stat(labels)
    sims = np.empty(N_PERM)
    idx = np.arange(len(units))
    for k in range(N_PERM):
        m = np.zeros(len(units), bool)
        m[rng.choice(idx, n1, replace=False)] = True
        sims[k] = stat(m)
    p = (np.sum(sims >= obs) + 1) / (N_PERM + 1)
    return obs, p, sims, len(units), n1


def temporal_stability(store, rng):
    rows = []
    for scheme in SCHEMES:
        pooled_tr = sum(store[(scheme, g, "treino", "todas")] for g in GENOS)
        for g in GENOS:
            other = GENOS[1] if g == GENOS[0] else GENOS[0]
            Ctr, Cte = store[(scheme, g, "treino", "todas")], store[(scheme, g, "teste", "todas")]
            ntr, nte = Ctr.sum(1), Cte.sum(1)
            R = [i for i in range(5) if ntr[i] > 0 and nte[i] > 0]
            Ptr = Ctr[R] / ntr[R, None]; Pte = Cte[R] / nte[R, None]
            Fr = float(np.linalg.norm(Ptr - Pte))
            sims = np.zeros(N_BOOT)
            for p_row, k in zip(Ptr, nte[R]):
                s = rng.multinomial(int(k), p_row, size=N_BOOT) / k
                sims += ((s - p_row) ** 2).sum(1)
            sims = np.sqrt(sims)
            p_boot = (np.sum(sims >= Fr) + 1) / (N_BOOT + 1)
            rows.append({
                "regra_estados": scheme, "experimento": g,
                "n_transicoes_treino": int(Ctr.sum()), "n_transicoes_teste": int(Cte.sum()),
                "linhas_comparaveis": len(R), "frobenius_linhas_comparaveis": Fr,
                "p_bootstrap_parametrico": p_boot,
                "logscore_matriz_propria": mean_logscore(smooth(Ctr), Cte),
                "logscore_matriz_outro_experimento": mean_logscore(smooth(store[(scheme, other, "treino", "todas")]), Cte),
                "logscore_matriz_agrupada": mean_logscore(smooth(pooled_tr), Cte),
                "logscore_uniforme": float(np.log(0.2)),
            })
    return pd.DataFrame(rows)


def chain_table(store):
    rows = []
    for scheme in SCHEMES:
        for g in GENOS:
            C = store[(scheme, g, "treino", "todas")]
            r = chain_props(C)
            for i in range(5):
                rows.append({
                    "regra_estados": scheme, "experimento": g, "estado": i + 1, "nome": STATE_NAMES[i + 1],
                    "cadeia_irredutivel": (r["irredutivel"] if r else False) if r else False,
                    "P_permanencia": (r["permanencia"][i] if r else np.nan),
                    "pi_estacionaria": (r["pi"][i] if r else np.nan),
                    "passos_medios_ate_estado5": (r["mfpt5"][i] if r else np.nan),
                    "nota": "" if r else "linha vazia: propriedades não definidas",
                })
    return pd.DataFrame(rows)


def cut_decomposition(trans):
    rows = []
    for scheme in SCHEMES:
        for g in GENOS:
            for s in SPLITS:
                t = trans[(trans.experimento == g) & (trans.amostra == s)]
                to1 = t[t[f"estado_{scheme}_para"] == 1]
                from5 = t[t[f"estado_{scheme}_de"] == 5]
                s5 = from5[from5[f"estado_{scheme}_para"] == 5]
                s5_cut = from5[from5.transicao_com_corte]
                rows.append({
                    "regra_estados": scheme, "experimento": g, "amostra": s,
                    "n_transicoes": len(t), "n_com_corte": int(t.transicao_com_corte.sum()),
                    "intervalo_medio_dias": t.intervalo_dias.mean(), "intervalo_mediano_dias": t.intervalo_dias.median(),
                    "intervalo_max_dias": t.intervalo_dias.max(),
                    "n_para_estado1": len(to1), "n_para_estado1_por_corte": int(to1.transicao_com_corte.sum()),
                    "n_saidas_estado5": len(from5), "n_saidas_estado5_por_corte": len(s5_cut),
                    "n_permanencias_estado5": len(s5),
                })
    return pd.DataFrame(rows)


def cycle_metrics(obs, mows, cwad):
    rows = []
    for key, d in obs.groupby("ciclo"):
        d = d.sort_values("data").reset_index(drop=True)
        g = d.experimento.iloc[0]
        y = int(d.year.iloc[0])
        cut = d.corte_ref.iloc[0]
        if len(d) < 3:
            rows.append({"ciclo": key, "experimento": g, "ano": y, "corte_inicio": cut, "n_obs": len(d),
                         "usado_nas_metricas": False})
            continue
        imax = int(d.LAI.values.argmax())
        imin = int(d.LAI.values[: imax + 1].argmin()) if imax > 0 else 0
        # último mínimo antes do máximo (início da rebrota)
        mn = d.LAI.values[: imax + 1].min()
        imin = int(np.where(d.LAI.values[: imax + 1] == mn)[0][-1])
        dt = (d.data.iloc[imax] - d.data.iloc[imin]).days
        # corte que ENCERRA o segmento (no mesmo ano); sem ele o ciclo é interrompido (dormência/fim da série)
        nxt = [c for c in mows[g] if c > cut and c.year == y]
        mes0 = d.data.iloc[imin].month
        estacao = 4 <= mes0 <= 10  # rebrota na estação principal de crescimento
        rebrota = bool(dt >= 7 and d.LAI.max() - mn >= 1.0 and estacao and nxt)
        row = {"ciclo": key, "experimento": g, "ano": y, "corte_inicio": cut, "n_obs": len(d),
               "mes_inicio_rebrota": mes0, "estacao_principal": estacao, "encerrado_por_corte": bool(nxt),
               "LAI_min_inicio": mn, "LAI_max": d.LAI.max(), "dias_min_ate_max": dt,
               "taxa_rebrota_LAI_dia": (d.LAI.max() - mn) / dt if rebrota else np.nan,
               "usado_nas_metricas": rebrota}
        if rebrota:
            seg = d.iloc[imin: imax + 1]
            for thr in (2.5, 3.8, 5.0):
                tt = np.nan
                x, yv = seg.data.map(lambda z: z.toordinal()).values, seg.LAI.values
                for i in range(1, len(seg)):
                    if yv[i - 1] < thr <= yv[i]:
                        tt = (x[i - 1] - x[0]) + (thr - yv[i - 1]) / (yv[i] - yv[i - 1]) * (x[i] - x[i - 1])
                        break
                row[f"dias_ate_LAI_{thr}"] = tt
        # biomassa (CWAD) do corte que encerra este ciclo
        if nxt:
            c_end = nxt[0]
            cw = cwad[(cwad.experimento == g) & (cwad.data >= c_end - timedelta(days=3)) & (cwad.data <= c_end)]
            row["corte_fim"] = c_end
            row["CWAD_kg_ha_no_corte"] = float(cw.CWAD.iloc[0]) if len(cw) else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["experimento", "corte_inicio", "ano"]).reset_index(drop=True)


def cycle_summary(cyc):
    u = cyc[cyc.usado_nas_metricas == True]
    rows = []
    for var in ["taxa_rebrota_LAI_dia", "LAI_max", "dias_ate_LAI_2.5", "dias_ate_LAI_3.8", "dias_ate_LAI_5.0", "CWAD_kg_ha_no_corte"]:
        if var not in u:
            continue
        for g in GENOS:
            v = u.loc[u.experimento == g, var].dropna()
            rows.append({"variavel": var, "experimento": g, "n": len(v), "media": v.mean() if len(v) else np.nan,
                         "dp": v.std(ddof=1) if len(v) > 1 else np.nan, "mediana": v.median() if len(v) else np.nan})
        a = u.loc[u.experimento == GENOS[0], var].dropna(); b = u.loc[u.experimento == GENOS[1], var].dropna()
        if len(a) >= 3 and len(b) >= 3:
            p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        else:
            p = np.nan
        rows.append({"variavel": var, "experimento": "comparacao_Mann_Whitney", "n": len(a) + len(b),
                     "media": np.nan, "dp": np.nan, "mediana": np.nan, "p_valor": p})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- contexto (fontes pesquisadas)
CONTEXT = [
    # indicador, valor, unidade, ano, fonte, url, observação
    ("Área de alfafa e misturas com gramíneas colhida para feno — EUA", 14.68, "milhões de acres (≈5,94 Mha)", 2025,
     "Hay & Forage Grower, com dados USDA/NASS", "https://hayandforage.com/article-5581-Slow-and-steady-hay-gains-in-2025.html",
     "Inclui misturas alfafa-gramíneas."),
    ("Produção de feno de alfafa e misturas — EUA", 50.21, "milhões de toneladas (USDA usa short tons)", 2025,
     "Hay & Forage Grower, com dados USDA/NASS", "https://hayandforage.com/article-5581-Slow-and-steady-hay-gains-in-2025.html",
     "Produtividade 3,42 t/acre; alta de <1% sobre 2024."),
    ("Área de alfafa — Argentina", 4_000_000, "ha", 2020, "ANBA (cita área da América Latina)",
     "https://anba.com.br/en/alfalfa-in-demand-on-international-markets/", "Maior área da América Latina."),
    ("Área de alfafa — Chile", 170_000, "ha", 2020, "ANBA", "https://anba.com.br/en/alfalfa-in-demand-on-international-markets/", ""),
    ("Área de alfafa — Peru", 120_000, "ha", 2020, "ANBA", "https://anba.com.br/en/alfalfa-in-demand-on-international-markets/", ""),
    ("Área de alfafa — Uruguai", 70_000, "ha", 2020, "ANBA", "https://anba.com.br/en/alfalfa-in-demand-on-international-markets/", ""),
    ("Área de alfafa — Brasil", 40_000, "ha", 2020, "ANBA; Embrapa (Tipificação de feno de alfafa para o Brasil)",
     "https://anba.com.br/en/alfalfa-in-demand-on-international-markets/", "Concentrada no RS e PR; Embrapa diz 'cerca de 40 mil ha ou menos'."),
    ("Área de alfafa — Aragón (Espanha)", 57_000, "ha (mais de)", 2020, "SomosAlfalfa (matéria de jul/2020)",
     "https://somosalfalfa.com.ar/negocios-noticias/espana-afianzada-como-primer-fabricante-y-exportador-mundial-de-alfalfa-deshidratada/",
     "Principalmente no Vale do Ebro; ~750 mil t. A mesma matéria cita Cataluña: 23.205 ha, ~250 mil t."),
    ("Exportações espanholas de alfafa desidratada", 750_000, "t (mais de)", 2024, "Diario del Campo",
     "https://diariodelcampo.com/las-exportaciones-espanolas-de-alfalfa-aumentaron-un-quince-por-ciento-en-2024-con-arabia-emiratos-y-china-como-grandes-compradores/",
     "Alta de ~15%; cerca de 2/3 em fardos e o restante em pellets. Principais destinos: Arábia Saudita, Emirados Árabes e China."),
    ("Comércio mundial de feno de alfafa — valor", 2.7, "bilhões US$", 2016, "ANBA",
     "https://anba.com.br/en/alfalfa-in-demand-on-international-markets/", "Demanda de ~6 milhões de t de feno de alfafa. Embrapa cita ~6,5 Mt de alfafa dentro de ~8,5 Mt de feno comercializado."),
    ("Participação de 5 importadores no consumo mundial de feno de alfafa", 84, "%", 2020, "ANBA",
     "https://anba.com.br/en/alfalfa-in-demand-on-international-markets/", "Arábia Saudita, Emirados, China, Japão e Coreia do Sul."),
    ("Área somada de EUA, UE, Argentina, Rússia, Canadá e Austrália", 32, "milhões de ha", 2020, "ANBA",
     "https://anba.com.br/en/alfalfa-in-demand-on-international-markets/", ""),
    ("Preço internacional do feno de alfafa (qualidade 'Good' a 'Supreme')", "220–280", "US$/t", 2023,
     "Embrapa (Tipificação de feno de alfafa para o Brasil)",
     "https://www.infoteca.cnptia.embrapa.br/infoteca/bitstream/doc/1165879/1/Tipificacao-de-feno-de-alfafa-para-o-Brasil.pdf", "Cotação de julho/2023."),
    ("Potencial produtivo com manejo e irrigação", 25, "t MS/ha/ano", np.nan, "Embrapa (Tipificação de feno de alfafa para o Brasil)",
     "https://www.infoteca.cnptia.embrapa.br/infoteca/bitstream/doc/1165879/1/Tipificacao-de-feno-de-alfafa-para-o-Brasil.pdf", "≈30 t/ha/ano de feno."),
    ("Proteína bruta da alfafa", "15–22", "%", np.nan, "NAAIC (Importance of Alfalfa)", "https://naaic.org/Alfalfa/Importance.html", "Uso principal: vacas leiteiras."),
    ("Modelo CROPGRO-Alfafa no Vale do Ebro — colheitas avaliadas", 58, "colheitas (990–4.617 kg/ha)", 2018,
     "Malik et al., Agronomy Journal 110(5):1777-1790", "https://doi.org/10.2134/agronj2017.12.0680",
     "LAI observado entre 0,1 e 6,7; RMSE 760 kg/ha; d=0,75 (produtividade) e 0,71 (LAI); 6 parcelas irrigadas por aspersão."),
]


def context_table() -> pd.DataFrame:
    return pd.DataFrame(CONTEXT, columns=["indicador", "valor", "unidade", "ano", "fonte", "url", "observacao"])


# ----------------------------------------------------------------------------- figuras
def make_figures(obs, mows, trans, store, cyc, perm_sims, perm_obs, stab, ctx):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 110, "savefig.dpi": 170, "savefig.bbox": "tight"})
    COL = {"AGZG1219": "#1b7f5a", "AGZG1501": "#c2570c"}
    BAND = ["#f4f1ea", "#eaf3ec", "#dbece0", "#c9e3d1", "#b6d9c1"]

    # 01 — série temporal com cortes
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.4), sharey=True)
    for ax, g in zip(axes, GENOS):
        d = obs[obs.experimento == g].sort_values("data")
        lo = [0] + EDGES
        hi = EDGES + [7.2]
        for l, h, c in zip(lo, hi, BAND):
            ax.axhspan(l, h, color=c, alpha=0.7, lw=0)
        for yr in (2016, 2017):
            dy = d[d.year == yr]
            ax.plot(dy.data, dy.LAI, "-o", ms=4, lw=1.4, color=COL[g])
        for c in mows[g]:
            if c.year in (2016, 2017):
                ax.axvline(c, color="#444", lw=0.8, ls="--")
        ax.set_title(f"{g} — LAI observado (linhas tracejadas = cortes)", loc="left", fontsize=10)
        ax.set_ylabel("LAI")
        ax.set_ylim(0, 7.2)
    fig.tight_layout(); fig.savefig(FIG / "01_serie_LAI_cortes.png"); plt.close(fig)

    # 02 — trajetórias por ciclo
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, g in zip(axes, GENOS):
        for key, d in obs[obs.experimento == g].groupby("ciclo"):
            d = d.sort_values("data")
            if len(d) < 3:
                continue
            x = (d.data - d.data.iloc[0]).map(lambda z: z.days)
            ax.plot(x, d.LAI, "-o", ms=3, lw=1, alpha=0.8, color=COL[g])
        for e in EDGES:
            ax.axhline(e, color="#999", lw=0.6, ls=":")
        ax.set_title(g, loc="left"); ax.set_xlabel("dias desde a primeira observação do ciclo")
    axes[0].set_ylabel("LAI")
    fig.tight_layout(); fig.savefig(FIG / "02_trajetorias_por_ciclo.png"); plt.close(fig)

    # 03 — heatmaps das matrizes de treino
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.6))
    for r, scheme in enumerate(["csv", "resumo"]):
        for c, g in enumerate(GENOS):
            C = store[(scheme, g, "treino", "todas")]; P = to_prob(C)
            ax = axes[r, c]
            ax.imshow(P, cmap="Greens", vmin=0, vmax=1)
            for i in range(5):
                for j in range(5):
                    if C[i].sum() == 0:
                        txt = "–"
                    else:
                        txt = f"{P[i, j]:.2f}\n(n={C[i, j]})" if C[i, j] else "·"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=7.5,
                            color="white" if P[i, j] > 0.55 else "#222")
            ax.set_xticks(range(5)); ax.set_yticks(range(5))
            ax.set_xticklabels(range(1, 6)); ax.set_yticklabels(range(1, 6))
            ax.set_title(f"{g} — {'estado do CSV' if scheme == 'csv' else 'limiares do resumo'}", fontsize=9.5, loc="left")
            ax.set_xlabel("estado seguinte"); ax.set_ylabel("estado atual")
            ax.spines[:].set_visible(False)
    fig.tight_layout(); fig.savefig(FIG / "03_matrizes_treino_2016.png"); plt.close(fig)

    # 04 — floresta de IC (P12 e P55)
    mt = pd.read_csv(OUT / "matrizes_transicao_long.csv")
    cells = [("AGZG1219", 1, 2), ("AGZG1501", 1, 2), ("AGZG1219", 5, 5), ("AGZG1501", 5, 5)]
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    y = 0
    ticks, labs = [], []
    for g, i, j in cells:
        for scheme, off in (("csv", 0.15), ("resumo", -0.15)):
            r = mt[(mt.regra_estados == scheme) & (mt.experimento == g) & (mt.amostra == "treino") &
                   (mt.subconjunto == "todas") & (mt.de == i) & (mt.para == j)].iloc[0]
            ax.errorbar(r.prob, y + off, xerr=[[r.prob - r.ic95_inf], [r.ic95_sup - r.prob]], fmt="o" if scheme == "csv" else "s",
                        color=COL[g], capsize=3, ms=5, mfc=COL[g] if scheme == "csv" else "white")
            ax.text(1.02, y + off, f"{r.contagem}/{r.n_linha}", va="center", fontsize=8)
        ticks.append(y); labs.append(f"{g}\nP{i}{j}"); y += 1
    ax.set_yticks(ticks); ax.set_yticklabels(labs); ax.set_xlim(0, 1.12); ax.invert_yaxis()
    ax.set_xlabel("probabilidade de transição (IC 95% de Clopper-Pearson)")
    ax.set_title("Os intervalos dos dois experimentos se sobrepõem  (● estado do CSV, □ limiares do resumo)", loc="left", fontsize=9.5)
    fig.tight_layout(); fig.savefig(FIG / "04_incerteza_P12_P55.png"); plt.close(fig)

    # 05 — nulo de permutação
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, scheme in zip(axes, ["csv", "resumo"]):
        sims = perm_sims[scheme]
        ax.hist(sims, bins=40, color="#9fbfae")
        ax.axvline(perm_obs[scheme][0], color="#c2570c", lw=2)
        ax.set_title(f"{'CSV' if scheme == 'csv' else 'Limiares do resumo'}: F obs = {perm_obs[scheme][0]:.3f}, p = {perm_obs[scheme][1]:.3f}",
                     fontsize=9.5, loc="left")
        ax.set_xlabel("distância de Frobenius sob permutação de ciclos")
    fig.tight_layout(); fig.savefig(FIG / "05_permutacao_frobenius.png"); plt.close(fig)

    # 06 — escore preditivo
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, scheme in zip(axes, ["csv", "resumo"]):
        s = stab[stab.regra_estados == scheme]
        x = np.arange(len(s)); w = 0.2
        series = [("logscore_matriz_propria", "matriz própria (2016)"), ("logscore_matriz_outro_experimento", "matriz do outro experimento"),
                  ("logscore_matriz_agrupada", "matriz agrupada"), ("logscore_uniforme", "uniforme (1/5)")]
        for k, (col, lab) in enumerate(series):
            ax.bar(x + (k - 1.5) * w, s[col], w, label=lab, color=["#1b7f5a", "#c2570c", "#7a8cc2", "#bbb"][k])
        ax.set_xticks(x); ax.set_xticklabels(s.experimento)
        ax.set_title(f"{'CSV' if scheme == 'csv' else 'Limiares do resumo'} — log-verossimilhança média em 2017", fontsize=9.5, loc="left")
    axes[0].set_ylabel("log-score por transição (maior = melhor)")
    axes[0].legend(fontsize=7.5, loc="lower left")
    fig.tight_layout(); fig.savefig(FIG / "06_validacao_preditiva_2017.png"); plt.close(fig)

    # 07 — contexto: áreas
    area = {
        "EUA (alfafa+mistura, 2025)": 14.68 * 0.404686 * 1e6,
        "Argentina (2020)": 4_000_000, "Chile (2020)": 170_000, "Peru (2020)": 120_000,
        "Uruguai (2020)": 70_000, "Aragón-ES (2020)": 57_000, "Brasil (2020)": 40_000,
    }
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    names = list(area); vals = [area[k] for k in names]
    colors = ["#c2570c" if "Aragón" in n else "#1b7f5a" if "Brasil" in n else "#8aa899" for n in names]
    ax.barh(names[::-1], vals[::-1], color=colors[::-1]); ax.set_xscale("log")
    for n, v in zip(names[::-1], vals[::-1]):
        ax.text(v * 1.08, n, f"{v/1e3:,.0f} mil ha".replace(",", "."), va="center", fontsize=8)
    ax.set_xlim(2e4, 3e7); ax.set_xlabel("área cultivada (ha, escala log)")
    ax.set_title("Área de alfafa: Brasil e Aragón (local dos dados) estão na mesma ordem de grandeza", loc="left", fontsize=9.2)
    fig.text(0.01, -0.03, "Fontes: Hay & Forage (USDA), ANBA 2020, SomosAlfalfa 2020. Anos e metodologias diferem — comparação de ordem de grandeza.", fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "07_contexto_area_alfafa.png"); plt.close(fig)


# ----------------------------------------------------------------------------- principal
def md_table(df: pd.DataFrame, digits=4) -> str:
    d = df.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.{digits}f}")
    head = "| " + " | ".join(map(str, d.columns)) + " |\n| " + " | ".join("---" for _ in d.columns) + " |\n"
    return head + "\n".join("| " + " | ".join(map(str, r)) + " |" for r in d.values)


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)

    obs, mows, cwad = build_observations()
    obs, chk = check_against_csv(obs)
    trans = build_transitions(obs)
    mt, store = matrices_tables(trans)
    fro_rows = frobenius_tables(store, rng)

    perm_sims, perm_obs = {}, {}
    perm_rows = []
    for scheme in SCHEMES:
        o, p, sims, nu, n1 = permutation_test(trans, scheme, rng)
        perm_sims[scheme], perm_obs[scheme] = sims, (o, p)
        perm_rows.append({"regra_estados": scheme, "frobenius_observado": o, "p_permutacao_por_ciclo": p,
                          "n_ciclos_treino": nu, "n_ciclos_AGZG1219": n1, "n_ciclos_AGZG1501": nu - n1, "n_permutacoes": N_PERM})
        for r in fro_rows:
            if r["tipo"] == "entre_experimentos" and r["regra_estados"] == scheme and "(treino)" in r["comparacao"]:
                r["p_permutacao"] = p
                r["nota"] = "p por permutação de ciclos de rebrota"
    stab = temporal_stability(store, rng)
    for r in fro_rows:
        if r["tipo"] == "treino_vs_teste":
            g = r["comparacao"].split(":")[0]
            s = stab[(stab.regra_estados == r["regra_estados"]) & (stab.experimento == g)].iloc[0]
            r["p_permutacao"] = s["p_bootstrap_parametrico"]
            r["nota"] = f"p = bootstrap paramétrico em {int(s['linhas_comparaveis'])} linhas comparáveis (F restrito = {s['frobenius_linhas_comparaveis']:.3f})"
    chain = chain_table(store)
    cutdec = cut_decomposition(trans)
    cyc = cycle_metrics(obs, mows, cwad)
    csum = cycle_summary(cyc)
    ctx = context_table()

    # sensibilidade dos números do resumo
    sens = []
    for scheme in SCHEMES:
        for g, (i, j) in [("AGZG1219", (1, 2)), ("AGZG1501", (5, 5))]:
            r = mt[(mt.regra_estados == scheme) & (mt.experimento == g) & (mt.amostra == "treino") &
                   (mt.subconjunto == "todas") & (mt.de == i) & (mt.para == j)].iloc[0]
            sens.append({"regra_estados": scheme, "metrica": f"P{i}{j} {g} (treino)", "valor": r.prob,
                         "contagem": f"{r.contagem}/{r.n_linha}", "ic95_inf": r.ic95_inf, "ic95_sup": r.ic95_sup})
        F = frob(to_prob(store[(scheme, GENOS[0], "treino", "todas")]), to_prob(store[(scheme, GENOS[1], "treino", "todas")]))
        sens.append({"regra_estados": scheme, "metrica": "Frobenius entre experimentos (treino)", "valor": F,
                     "contagem": "", "ic95_inf": np.nan, "ic95_sup": np.nan})
    sens = pd.DataFrame(sens)

    # gravação
    obs_out = obs.drop(columns=["DATE"]).copy()
    obs_out.to_csv(OUT / "observacoes_lai.csv", index=False)
    chk.to_csv(OUT / "verificacao_csv_vs_bruto.csv", index=False)
    trans.to_csv(OUT / "transicoes.csv", index=False)
    mt.to_csv(OUT / "matrizes_transicao_long.csv", index=False)
    pd.DataFrame(fro_rows).to_csv(OUT / "distancias_frobenius.csv", index=False)
    pd.DataFrame(perm_rows).to_csv(OUT / "teste_permutacao_entre_experimentos.csv", index=False)
    stab.to_csv(OUT / "estabilidade_temporal_2017.csv", index=False)
    chain.to_csv(OUT / "propriedades_cadeia.csv", index=False)
    cutdec.to_csv(OUT / "decomposicao_cortes.csv", index=False)
    cyc.to_csv(OUT / "ciclos_rebrota.csv", index=False)
    csum.to_csv(OUT / "resumo_ciclos.csv", index=False)
    sens.to_csv(OUT / "sensibilidade_regra_de_estados.csv", index=False)
    ctx.to_csv(OUT / "contexto_alfafa_fontes.csv", index=False)
    pd.DataFrame([{"estado": k, "nome": v, "regra_resumo": r} for (k, v), r in
                  zip(STATE_NAMES.items(), ["LAI < 1,0", "1,0 ≤ LAI < 2,5", "2,5 ≤ LAI < 3,8", "3,8 ≤ LAI < 5,0", "LAI ≥ 5,0"])]
                 ).to_csv(OUT / "definicao_estados.csv", index=False)

    make_figures(obs, mows, trans, store, cyc, perm_sims, perm_obs, stab, ctx)

    # snapshot
    key = {
        "P12_AGZG1219_csv": float(sens[(sens.regra_estados == "csv") & sens.metrica.str.startswith("P12")].valor.iloc[0]),
        "P55_AGZG1501_csv": float(sens[(sens.regra_estados == "csv") & sens.metrica.str.startswith("P55")].valor.iloc[0]),
        "Frobenius_treino_csv": float(sens[(sens.regra_estados == "csv") & sens.metrica.str.startswith("Frobenius")].valor.iloc[0]),
        "P12_AGZG1219_resumo": float(sens[(sens.regra_estados == "resumo") & sens.metrica.str.startswith("P12")].valor.iloc[0]),
        "P55_AGZG1501_resumo": float(sens[(sens.regra_estados == "resumo") & sens.metrica.str.startswith("P55")].valor.iloc[0]),
        "Frobenius_treino_resumo": float(sens[(sens.regra_estados == "resumo") & sens.metrica.str.startswith("Frobenius")].valor.iloc[0]),
        "p_permutacao_csv": perm_obs["csv"][1], "p_permutacao_resumo": perm_obs["resumo"][1],
        "n_observacoes": len(obs),
    }
    snap = {
        "seed": SEED, "n_permutacoes": N_PERM, "n_bootstrap": N_BOOT, "limiares_LAI": EDGES,
        "genotipos": GENOS, "divisao": {"treino": 2016, "teste": 2017},
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "arquivos_brutos_sha256": {p.name: sha256(p) for p in sorted(RAW.iterdir())},
        "metricas_chave": key,
    }
    (ROOT / "config_snapshot.json").write_text(json.dumps(snap, indent=2, ensure_ascii=False))

    # NotebookLM
    NLM.mkdir(exist_ok=True)
    for f in NLM.glob("*.csv"):
        f.unlink()
    order = ["verificacao_csv_vs_bruto", "definicao_estados", "observacoes_lai", "transicoes", "matrizes_transicao_long",
             "sensibilidade_regra_de_estados", "distancias_frobenius", "teste_permutacao_entre_experimentos",
             "estabilidade_temporal_2017", "propriedades_cadeia", "decomposicao_cortes", "ciclos_rebrota", "resumo_ciclos",
             "contexto_alfafa_fontes"]
    for k, n in enumerate(order, start=1):
        shutil.copy(OUT / f"{n}.csv", NLM / f"{k:02d}_{n}.csv")

    # matrizes de treino em texto para o master
    mats = []
    for scheme in ["csv", "resumo"]:
        for g in GENOS:
            C = store[(scheme, g, "treino", "todas")]
            P = pd.DataFrame(to_prob(C), index=[f"E{i}" for i in range(1, 6)], columns=[f"E{i}" for i in range(1, 6)]).round(3)
            Cn = pd.DataFrame(C, index=P.index, columns=P.columns)
            mats.append(f"### {g} — treino 2016 — regra: {SCHEME_LABEL[scheme]}\n\nProbabilidades:\n\n{md_table(P.reset_index().rename(columns={'index': 'de/para'}), 3)}\n\nContagens:\n\n{md_table(Cn.reset_index().rename(columns={'index': 'de/para'}), 0)}\n")
    master = f"""# MASTER DE DADOS PARA NOTEBOOKLM — Modelagem markoviana do LAI em alfafa (IX CBRG)

Repositório: https://github.com/DaviHrAmancio/GeneticaAlfafa

## Regras de integridade
Classificar cada conteúdo como: DADO EMPÍRICO (LAI observado, cortes) · PARÂMETRO ESTIMADO (matrizes, IC) ·
TESTE ESTATÍSTICO (permutação, bootstrap) · CONTEXTO DE LITERATURA (fontes públicas citadas). Não inventar números.
Os experimentos AGZG1219 e AGZG1501 são códigos de experimento DSSAT; confirmar com o orientador se representam genótipos distintos.

## Dado de origem
{len(obs)} observações de LAI extraídas dos arquivos .ALT (campo LAID); cortes dos arquivos .MOW; experimentos na Espanha (Aragón/Vale do Ebro).
Treino = 2016; teste = 2017. Verificação do CSV contra o bruto:

{md_table(chk, 4)}

## Duas regras de discretização (importante)
- `csv`: coluna `estado` do CSV preparado pelo ChatGPT. Reproduz exatamente os números do resumo submetido.
- `resumo`: limiares do resumo (1,0 / 2,5 / 3,8 / 5,0) aplicados ao LAI bruto. {int(obs.estado_difere.sum())} das {len(obs)} observações mudam de estado.

## Números do resumo e sensibilidade

{md_table(sens, 4)}

## Teste de diferença entre experimentos (permutação por ciclo de rebrota, treino)

{md_table(pd.DataFrame(perm_rows), 4)}

## Estabilidade temporal 2016→2017 e capacidade preditiva

{md_table(stab, 4)}

## Decomposição por cortes (a dinâmica é dominada pelo manejo de corte)

{md_table(cutdec, 2)}

## Métricas de rebrota por ciclo (sem discretização)

{md_table(csum, 4)}

## Matrizes de transição (treino)

{chr(10).join(mats)}

## Contexto da alfafa (fontes públicas)

{md_table(ctx[['indicador','valor','unidade','ano','fonte']], 4)}
"""
    (NLM / "00_NOTEBOOKLM_MASTER.md").write_text(master, encoding="utf-8")
    for d in ["METODOLOGIA", "TRANSPARENCIA_E_LIMITACOES", "CONTEXTO_ALFAFA", "MANIFESTO_DE_EVIDENCIAS"]:
        src = ROOT / "docs" / f"{d}.md"
        if src.exists():
            shutil.copy(src, NLM / f"{15 + ['METODOLOGIA','TRANSPARENCIA_E_LIMITACOES','CONTEXTO_ALFAFA','MANIFESTO_DE_EVIDENCIAS'].index(d)}_{d}.md")

    print("OK — chaves:")
    print(json.dumps(key, indent=2))


if __name__ == "__main__":
    main()
