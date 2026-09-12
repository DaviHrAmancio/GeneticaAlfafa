"""
Explorador Didático — Dinâmica Fenotípica da Alfafa (Medicago sativa L.)
===========================================================================
Projeto: IX CBRG 2026 — Davi Henrique de Souza Amâncio & Prof. Leonardo Talavera Campos
IFC Campus Camboriú

Este programa reúne DUAS fontes de dados em uma interface só, organizada em abas:

  Aba 1 — "Sua Análise (Markov)"
      A análise original do artigo: LAI de AGZG1219 e AGZG1501, discretizado em
      5 estados fenotípicos, com a matriz de transição de Markov (dados do seu
      arquivo alfafa_markov_cbrg.db).

  Aba 2 — "Base do Professor (DSSAT)"
      Os 5 experimentos originais do pacote DSSAT/dssat-csm-data (Alfalfa.zip):
      AGZG1219, AGZG1229, AGZG1501, AGZG1502 (Espanha) e ONCF1401 (Ottawa,
      Canadá). Nem todos têm LAI medido — alguns só têm biomassa (CWAD/HERB).

  Aba 3 — "Comparativo Multi-Local"
      Sobreposição da biomassa acumulada (CWAD, a única variável presente em
      TODOS os experimentos) para comparar visualmente o efeito do clima e do
      manejo entre Espanha e Canadá.

  Aba 5 — "Como o Markov Funciona"
      O "motor" do modelo, explicado por dentro, em 3 painéis:
        A) Construção da matriz: anda observação a observação pela série REAL de
           treino, contando cada transição ao vivo até formar a matriz final —
           e confere que bate com a matriz oficial do artigo.
        B) Diagrama de estados: grafo com os 5 estados e setas cuja espessura
           representa a probabilidade de transição.
        C) Simulação de Monte Carlo: "roda o dado" segundo a matriz, gera uma
           trajetória sintética passo a passo (com a régua de sorteio visível),
           e compara a frequência simulada de cada estado com a distribuição
           estacionária teórica (o que acontece no longuíssimo prazo).

  Aba 4 — "Metodologia e Alertas"
      Explicação botânica + observações importantes sobre os dados, para você
      levar em conta antes de apresentar/submeter ao IX CBRG.

COMO EXECUTAR
--------------
1. Deixe este arquivo, o "alfafa_markov_cbrg.db" e o "Alfalfa.zip"
   (ou a pasta "Alfalfa" já descompactada) na MESMA pasta.
2. Instale as dependências:
       pip install customtkinter matplotlib pandas numpy
3. Rode:
       python explorador_alfafa_v3.py
"""

import math
import random
import re
import sys
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

try:
    import customtkinter as ctk
except ImportError:
    sys.exit(
        "\nA biblioteca 'customtkinter' não está instalada.\n"
        "Instale com:  pip install customtkinter matplotlib pandas\n"
    )

import sqlite3

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "alfafa_markov_cbrg.db"
ZIP_PATH = BASE_DIR / "Alfalfa.zip"
ALFALFA_DIR = BASE_DIR / "Alfalfa"

DSSAT_EXPERIMENTS = ["AGZG1219", "AGZG1229", "AGZG1501", "AGZG1502", "ONCF1401"]

STATE_COLORS = {1: "#c0392b", 2: "#e67e22", 3: "#f1c40f", 4: "#2ecc71", 5: "#27ae60"}
STATE_NAMES = {
    1: "Pós-corte", 2: "Rebrota inicial", 3: "Crescimento intermediário",
    4: "Fechamento de dossel", 5: "Cobertura plena",
}
STATE_MECHANISM = {
    1: ("A parte aérea fotossinteticamente ativa foi quase totalmente removida "
        "pelo corte. A planta deixa de depender da fotossíntese corrente e passa "
        "a rebrotar às custas das RESERVAS de carboidratos e nitrogênio armazenadas "
        "na coroa (crown) e nas raízes."),
    2: ("Gemas dormentes da coroa emitem novas hastes com folhas trifolioladas "
        "pequenas. A fotossíntese recomeça, mas a planta ainda está 'pagando' "
        "energeticamente a reconstrução do aparato foliar."),
    3: ("Fase de maior taxa de crescimento relativo: a área foliar já autossustenta "
        "a planta energeticamente e a expansão vegetativa acelera."),
    4: ("As folhas de hastes/plantas vizinhas começam a se sobrepor (fechamento de "
        "dossel). O autossombreamento das folhas mais baixas já reduz eficiência."),
    5: ("LAI próximo do platô fisiológico. O ganho marginal de biomassa por folha "
        "adicional cai bastante — na prática, é por volta daqui que se corta."),
}

MOW_DATES_DB = {  # já usados na aba 1 (janela 2016-2017 do artigo)
    "AGZG1219": ["2016-04-27", "2016-06-07", "2016-07-22", "2016-08-31",
                 "2016-10-21", "2017-04-05"],
    "AGZG1501": ["2016-04-25", "2016-05-25", "2016-07-06", "2016-08-02",
                 "2016-09-06", "2016-11-07", "2017-04-12", "2017-05-23"],
}

VAR_LABELS = {
    "LAID": "LAI (m²/m²)", "CWAD": "Biomassa aérea seca (kg/ha)",
    "HERB": "Herbagem colhida no corte (kg/ha)", "CPROT": "Proteína bruta (%)",
    "VN%D": "N foliar (%)", "CHTD": "Altura do dossel (m)",
}

EXP_LABELS = {
    "AGZG1219": "AGZG1219 — Finca, Espanha (multi-talhão, 2012-2017)",
    "AGZG1229": "AGZG1229 — Cavero, Espanha (2012-2014, sem LAI)",
    "AGZG1501": "AGZG1501 — Montañana/Aragón, Espanha (irrigado, 2016-2017)",
    "AGZG1502": "AGZG1502 — Montañana/Aragón (mesmos dados de AGZG1501, outra regra de corte)",
    "ONCF1401": "ONCF1401 — Ottawa, Canadá (clima temperado frio, sem LAI)",
}

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("green")


# --------------------------------------------------------------------------------
# PARSER GENÉRICO DE ARQUIVOS DSSAT (.ALT / .MOW / .ALX)
# --------------------------------------------------------------------------------

def julian_to_date(token):
    token = str(int(float(token)))
    if len(token) <= 5:
        token = token.zfill(5)
        yy, ddd = int(token[:2]), int(token[2:])
        year = 2000 + yy if yy < 50 else 1900 + yy
    else:
        token = token.zfill(7)
        year, ddd = int(token[:4]), int(token[4:])
    return datetime(year, 1, 1) + timedelta(days=ddd - 1)


def parse_alt(path):
    lines = Path(path).read_text(errors="ignore").splitlines()
    records, header = [], None
    for line in lines:
        if line.strip().startswith("@TRNO"):
            header = line.replace("@", " ").split()
            continue
        if header and line.strip() and not line.strip().startswith(("!", "*")):
            parts = line.split()
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            try:
                dt = julian_to_date(row["DATE"])
            except Exception:
                continue
            for var in header:
                if var in ("TRNO", "DATE"):
                    continue
                try:
                    val = float(row[var])
                except ValueError:
                    continue
                if val == -99:
                    continue
                records.append({"date": dt, "variable": var, "value": val})
    df = pd.DataFrame.from_records(records)
    if not df.empty:
        df = df.sort_values("date").reset_index(drop=True)
    return df


def parse_mow(path):
    lines = Path(path).read_text(errors="ignore").splitlines()
    dates, header = [], None
    for line in lines:
        if line.strip().startswith("@TRNO"):
            header = line.split()
            continue
        if header and line.strip() and not line.strip().startswith("!"):
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                dates.append(julian_to_date(parts[1]))
            except Exception:
                pass
    return sorted(set(dates))


def parse_alx_meta(path):
    text = Path(path).read_text(errors="ignore")
    meta = {}
    m = re.search(r"\*EXP\.DETAILS:\s*(.+)", text)
    meta["titulo"] = m.group(1).strip() if m else path.stem
    m = re.search(r"@SITE\s*\n(.+)", text)
    meta["site"] = m.group(1).strip() if m else "?"
    meta["n_irrig_eventos"] = len(re.findall(r"\bIR00\d\b", text))
    n_kg = re.search(r"(\d+\s*kg\s*N/ha)", text)
    meta["adubacao"] = n_kg.group(1) if n_kg else "não informado no cabeçalho"
    return meta


def ensure_dssat_extracted():
    if ALFALFA_DIR.exists():
        return True
    if ZIP_PATH.exists():
        with zipfile.ZipFile(ZIP_PATH) as z:
            z.extractall(BASE_DIR)
        return ALFALFA_DIR.exists()
    return False


def load_dssat_experiments():
    data = {}
    if not ensure_dssat_extracted():
        return data
    for exp in DSSAT_EXPERIMENTS:
        alt, mow, alx = (ALFALFA_DIR / f"{exp}.ALT", ALFALFA_DIR / f"{exp}.MOW",
                         ALFALFA_DIR / f"{exp}.ALX")
        if not alt.exists():
            continue
        df = parse_alt(alt)
        mows = parse_mow(mow) if mow.exists() else []
        meta = parse_alx_meta(alx) if alx.exists() else {}
        data[exp] = {"df": df, "mow": mows, "meta": meta}
    return data


def load_markov_db():
    if not DB_PATH.exists():
        return None, None, None
    conn = sqlite3.connect(DB_PATH)
    obs = pd.read_sql_query(
        "SELECT data, year, doy, LAI, experimento, amostra, "
        "dias_desde_primeira_obs, estado FROM observacoes_lai ORDER BY experimento, data",
        conn,
    )
    trans = pd.read_sql_query("SELECT * FROM matrizes_transicao WHERE amostra='treino'", conn)
    conn.close()
    obs["data"] = pd.to_datetime(obs["data"])
    return obs, trans, True


# --------------------------------------------------------------------------------
# DIAGRAMA DE ESTADOS DE MARKOV (desenhado manualmente, sem depender de networkx)
# --------------------------------------------------------------------------------

def state_positions(radius=1.3):
    """5 estados dispostos em pentágono, na ordem 1→2→3→4→5→(corte)→1."""
    pos = {}
    for s in range(1, 6):
        ang = math.pi / 2 - (s - 1) * (2 * math.pi / 5)   # começa no topo, sentido horário
        pos[s] = (radius * math.cos(ang), radius * math.sin(ang))
    return pos


def draw_markov_diagram(ax, matrix_rows, highlight=None, min_prob=0.03, title=None):
    """matrix_rows: lista de dicts/linhas com estado_origem, estado_destino, probabilidade."""
    ax.clear()
    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(-2.1, 2.1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor("#1f1f1f")
    pos = state_positions()

    for row in matrix_rows:
        i, j, p = int(row["estado_origem"]), int(row["estado_destino"]), float(row["probabilidade"])
        if p < min_prob:
            continue
        xi, yi = pos[i]
        if i == j:
            dx, dy = xi / 1.3, yi / 1.3
            p1 = (xi + 0.18 * dx - 0.16 * dy, yi + 0.18 * dy + 0.16 * dx)
            p2 = (xi + 0.18 * dx + 0.16 * dy, yi + 0.18 * dy - 0.16 * dx)
            loop = mpatches.FancyArrowPatch(
                p1, p2, connectionstyle="arc3,rad=1.8", arrowstyle="-|>",
                mutation_scale=11, lw=0.6 + 7 * p, alpha=0.85, color="#e67e22", zorder=3,
            )
            ax.add_patch(loop)
            ax.text(xi + 0.55 * dx, yi + 0.55 * dy, f"{p*100:.0f}%", fontsize=8,
                     color="#f5c16c", ha="center", zorder=4)
        else:
            xj, yj = pos[j]
            rad = 0.22 if i < j else -0.30
            arrow = mpatches.FancyArrowPatch(
                (xi, yi), (xj, yj), connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                mutation_scale=13, lw=0.6 + 7 * p, alpha=0.80, color="#5dade2",
                shrinkA=20, shrinkB=20, zorder=2,
            )
            ax.add_patch(arrow)
            mx, my = (xi + xj) / 2, (yi + yj) / 2
            ox, oy = -(yj - yi), (xj - xi)
            norm = math.hypot(ox, oy) or 1
            mx += 0.18 * rad / abs(rad) * ox / norm if rad else 0
            my += 0.18 * rad / abs(rad) * oy / norm if rad else 0
            ax.text(mx, my, f"{p*100:.0f}%", fontsize=8, color="white", ha="center", zorder=4)

    for s, (x, y) in pos.items():
        edge = "white" if highlight == s else "#2a2a2a"
        lw = 3.2 if highlight == s else 1.2
        circ = mpatches.Circle((x, y), 0.32, facecolor=STATE_COLORS[s], edgecolor=edge, linewidth=lw, zorder=5)
        ax.add_patch(circ)
        ax.text(x, y, str(s), ha="center", va="center", fontsize=15, fontweight="bold",
                 color="#1a1a1a", zorder=6)
        ax.text(x * 1.55, y * 1.55, STATE_NAMES[s], ha="center", va="center", fontsize=7.5,
                 color="#cfd8dc", zorder=6)
    if title:
        ax.set_title(title, color="white", fontsize=11)


def stationary_distribution(P, iters=500):
    """Distribuição estacionária via iteração de potência (P deve ser 5x5, linhas somando 1)."""
    n = P.shape[0]
    v = np.ones(n) / n
    for _ in range(iters):
        v = v @ P
        s = v.sum()
        if s > 0:
            v = v / s
    return v


def matrix_df_to_numpy(sub_trans):
    """Converte linhas (estado_origem, estado_destino, probabilidade) num array 5x5,
    completando estados nunca observados como origem com auto-loop (permanece)."""
    P = np.zeros((5, 5))
    for _, r in sub_trans.iterrows():
        P[int(r["estado_origem"]) - 1, int(r["estado_destino"]) - 1] = r["probabilidade"]
    for i in range(5):
        if P[i].sum() == 0:
            P[i, i] = 1.0
    return P


# --------------------------------------------------------------------------------
# APLICATIVO PRINCIPAL
# --------------------------------------------------------------------------------

class ExploradorAlfafa(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Explorador Didático — Alfafa (IX CBRG)")
        self.geometry("1260x800")
        self.minsize(1040, 680)

        self.obs, self.trans, ok = load_markov_db()
        self.dssat = load_dssat_experiments()

        tabview = ctk.CTkTabview(self, width=1220, height=770)
        tabview.pack(fill="both", expand=True, padx=10, pady=10)

        t1 = tabview.add("① Sua Análise (Markov)")
        t2 = tabview.add("② Base do Professor (DSSAT)")
        t3 = tabview.add("③ Comparativo Multi-Local")
        t5 = tabview.add("⑤ Como o Markov Funciona")
        t4 = tabview.add("④ Metodologia e Alertas")

        if self.obs is not None:
            self._build_markov_tab(t1)
            self._build_howmarkov_tab(t5)
        else:
            ctk.CTkLabel(t1, text="alfafa_markov_cbrg.db não encontrado nesta pasta.").pack(pady=40)
            ctk.CTkLabel(t5, text="alfafa_markov_cbrg.db não encontrado nesta pasta.").pack(pady=40)

        if self.dssat:
            self._build_dssat_tab(t2)
            self._build_compare_tab(t3)
        else:
            for t in (t2, t3):
                ctk.CTkLabel(
                    t, text="Alfalfa.zip não encontrado nesta pasta (ou pasta 'Alfalfa' extraída)."
                ).pack(pady=40)

        self._build_info_tab(t4)

    # ============================================================ ABA 1 =====
    def _build_markov_tab(self, parent):
        parent.grid_columnconfigure(0, weight=0)
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        side = ctk.CTkFrame(parent, width=330, corner_radius=0)
        side.grid(row=0, column=0, sticky="nswe")
        side.grid_propagate(False)

        ctk.CTkLabel(side, text="Genótipo do artigo", font=ctk.CTkFont(size=14, weight="bold")).pack(
            anchor="w", padx=16, pady=(16, 4)
        )
        self.m_geno = ctk.CTkSegmentedButton(
            side, values=["AGZG1219", "AGZG1501"], command=lambda v: self._m_set_geno(v)
        )
        self.m_geno.set("AGZG1219")
        self.m_geno.pack(fill="x", padx=16, pady=(0, 12))

        self.m_play_btn = ctk.CTkButton(side, text="▶ Reproduzir", command=self._m_toggle_play)
        self.m_play_btn.pack(fill="x", padx=16, pady=(0, 8))

        self.m_slider = ctk.CTkSlider(side, from_=0, to=1, number_of_steps=1,
                                       command=lambda v: self._m_set_index(int(round(float(v)))))
        self.m_slider.pack(fill="x", padx=16, pady=(4, 10))

        self.m_date_label = ctk.CTkLabel(side, text="", font=ctk.CTkFont(size=14, weight="bold"))
        self.m_date_label.pack(anchor="w", padx=16)
        self.m_lai_label = ctk.CTkLabel(side, text="")
        self.m_lai_label.pack(anchor="w", padx=16, pady=(0, 10))

        self.m_state_card = ctk.CTkFrame(side, corner_radius=10)
        self.m_state_card.pack(fill="x", padx=16, pady=(0, 10))
        self.m_state_title = ctk.CTkLabel(self.m_state_card, text="", font=ctk.CTkFont(size=15, weight="bold"))
        self.m_state_title.pack(anchor="w", padx=12, pady=(10, 4))
        self.m_state_desc = ctk.CTkLabel(self.m_state_card, text="", justify="left", wraplength=280,
                                          font=ctk.CTkFont(size=11))
        self.m_state_desc.pack(anchor="w", padx=12, pady=(0, 10))

        ctk.CTkLabel(side, text="Probabilidade de transição (estado atual → amanhã)",
                     justify="left", wraplength=290, font=ctk.CTkFont(size=11)).pack(
            fill="x", padx=16, pady=(2, 4)
        )
        self.m_trans_box = ctk.CTkTextbox(side, height=110, font=ctk.CTkFont(size=12))
        self.m_trans_box.pack(fill="x", padx=16)
        self.m_trans_box.configure(state="disabled")

        main = ctk.CTkFrame(parent, fg_color="transparent")
        main.grid(row=0, column=1, sticky="nswe", padx=10, pady=10)
        main.grid_rowconfigure(0, weight=1)
        main.grid_columnconfigure(0, weight=1)
        self.m_fig, self.m_ax = plt.subplots(figsize=(7.6, 5.6), dpi=100)
        self.m_fig.patch.set_facecolor("#1a1a1a")
        self.m_canvas = FigureCanvasTkAgg(self.m_fig, master=main)
        self.m_canvas.get_tk_widget().grid(row=0, column=0, sticky="nswe")

        self.m_playing = False
        self._m_refresh_series()
        self._m_set_index(0)

    def _m_refresh_series(self):
        geno = self.m_geno.get()
        self.m_series = self.obs[self.obs.experimento == geno].sort_values("data").reset_index(drop=True)
        self.m_mow = pd.to_datetime(MOW_DATES_DB.get(geno, []))
        n = max(len(self.m_series) - 1, 1)
        self.m_slider.configure(from_=0, to=n, number_of_steps=n)

    def _m_set_geno(self, value):
        self._m_refresh_series()
        self._m_set_index(0)

    def _m_toggle_play(self):
        self.m_playing = not self.m_playing
        self.m_play_btn.configure(text="⏸ Pausar" if self.m_playing else "▶ Reproduzir")
        if self.m_playing:
            self._m_auto()

    def _m_auto(self):
        if not self.m_playing:
            return
        nxt = self.m_idx + 1
        if nxt >= len(self.m_series):
            self.m_playing = False
            self.m_play_btn.configure(text="▶ Reproduzir")
            return
        self._m_set_index(nxt)
        self.after(450, self._m_auto)

    def _m_set_index(self, idx):
        idx = int(max(0, min(idx, len(self.m_series) - 1)))
        self.m_idx = idx
        self.m_slider.set(idx)
        row = self.m_series.iloc[idx]
        estado = int(row["estado"])
        self.m_date_label.configure(text=row["data"].strftime("%d/%m/%Y") + f" · {row['amostra']}")
        self.m_lai_label.configure(text=f"LAI = {row['LAI']:.2f} m²/m²")
        self.m_state_card.configure(border_width=2, border_color=STATE_COLORS[estado])
        self.m_state_title.configure(text=f"Estado {estado} — {STATE_NAMES[estado]}",
                                      text_color=STATE_COLORS[estado])
        self.m_state_desc.configure(text=STATE_MECHANISM[estado])

        sub = self.trans[(self.trans.experimento == self.m_geno.get()) &
                          (self.trans.estado_origem == estado)].sort_values("estado_destino")
        self.m_trans_box.configure(state="normal")
        self.m_trans_box.delete("1.0", "end")
        for _, r in sub.iterrows():
            barra = "█" * int(round(r["probabilidade"] * 20))
            self.m_trans_box.insert("end", f"→ Est.{int(r['estado_destino'])}  {r['probabilidade']*100:5.1f}%  {barra}\n")
        self.m_trans_box.configure(state="disabled")

        self._m_draw()

    def _m_draw(self):
        ax = self.m_ax
        ax.clear()
        s = self.m_series
        ax.axhspan(0, 0.9, color=STATE_COLORS[1], alpha=0.10)
        ax.axhspan(0.9, 1.8, color=STATE_COLORS[2], alpha=0.10)
        ax.axhspan(1.8, 3.3, color=STATE_COLORS[3], alpha=0.10)
        ax.axhspan(3.3, 4.8, color=STATE_COLORS[4], alpha=0.10)
        ax.axhspan(4.8, max(s["LAI"].max() * 1.1, 6.5), color=STATE_COLORS[5], alpha=0.10)
        tr = s[s.amostra == "treino"]; te = s[s.amostra == "teste"]
        ax.plot(tr["data"], tr["LAI"], "-o", color="#5dade2", ms=4, lw=1.6, label="Treino (2016)")
        ax.plot(te["data"], te["LAI"], "--s", color="#f39c12", ms=5, lw=1.6, label="Teste (2017)")
        for md in self.m_mow:
            ax.axvline(md, color="#e74c3c", ls=":", lw=1, alpha=0.7)
        cur = s.iloc[self.m_idx]
        ax.scatter([cur["data"]], [cur["LAI"]], s=170, facecolor="none", edgecolor="white", lw=2.2, zorder=5)
        ax.set_facecolor("#1f1f1f")
        ax.set_title(f"{self.m_geno.get()} — LAI (linhas vermelhas = corte)", color="white", fontsize=11)
        ax.set_ylabel("LAI (m²/m²)", color="white")
        ax.tick_params(colors="white", labelrotation=25)
        for sp in ax.spines.values():
            sp.set_color("#555")
        ax.legend(fontsize=8, facecolor="#2a2a2a", labelcolor="white")
        self.m_fig.tight_layout()
        self.m_canvas.draw()

    # ============================================================ ABA 2 =====
    def _build_dssat_tab(self, parent):
        parent.grid_columnconfigure(0, weight=0)
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        side = ctk.CTkFrame(parent, width=340, corner_radius=0)
        side.grid(row=0, column=0, sticky="nswe")
        side.grid_propagate(False)

        ctk.CTkLabel(side, text="Experimento DSSAT", font=ctk.CTkFont(size=14, weight="bold")).pack(
            anchor="w", padx=16, pady=(16, 4)
        )
        self.d_exp_combo = ctk.CTkComboBox(
            side, values=[EXP_LABELS[e] for e in self.dssat.keys()],
            command=self._d_on_exp_change, width=300,
        )
        self.d_exp_combo.pack(padx=16, pady=(0, 12))

        self.d_meta_box = ctk.CTkTextbox(side, height=150, font=ctk.CTkFont(size=12))
        self.d_meta_box.pack(fill="x", padx=16, pady=(0, 12))
        self.d_meta_box.configure(state="disabled")

        ctk.CTkLabel(side, text="Variável a exibir").pack(anchor="w", padx=16)
        self.d_var_combo = ctk.CTkComboBox(side, values=[], command=lambda v: self._d_draw(), width=300)
        self.d_var_combo.pack(padx=16, pady=(2, 14))

        ctk.CTkLabel(
            side, text="Linhas vermelhas pontilhadas = eventos de corte (mow)\n"
            "registrados no arquivo .MOW original do experimento.",
            justify="left", wraplength=300, font=ctk.CTkFont(size=11), text_color="#b8c2cc",
        ).pack(fill="x", padx=16)

        main = ctk.CTkFrame(parent, fg_color="transparent")
        main.grid(row=0, column=1, sticky="nswe", padx=10, pady=10)
        main.grid_rowconfigure(0, weight=1)
        main.grid_columnconfigure(0, weight=1)
        self.d_fig, self.d_ax = plt.subplots(figsize=(7.6, 5.6), dpi=100)
        self.d_fig.patch.set_facecolor("#1a1a1a")
        self.d_canvas = FigureCanvasTkAgg(self.d_fig, master=main)
        self.d_canvas.get_tk_widget().grid(row=0, column=0, sticky="nswe")

        first = list(self.dssat.keys())[0]
        self.d_exp_combo.set(EXP_LABELS[first])
        self._d_on_exp_change(EXP_LABELS[first])

    def _d_current_exp_code(self):
        label = self.d_exp_combo.get()
        for code, lab in EXP_LABELS.items():
            if lab == label:
                return code
        return list(self.dssat.keys())[0]

    def _d_on_exp_change(self, _label):
        code = self._d_current_exp_code()
        info = self.dssat[code]
        df, mow, meta = info["df"], info["mow"], info["meta"]

        self.d_meta_box.configure(state="normal")
        self.d_meta_box.delete("1.0", "end")
        self.d_meta_box.insert("end", f"{meta.get('titulo','')}\n\n")
        self.d_meta_box.insert("end", f"Local: {meta.get('site','?')}\n")
        self.d_meta_box.insert("end", f"Eventos de irrigação no arquivo: {meta.get('n_irrig_eventos',0)}\n")
        self.d_meta_box.insert("end", f"Adubação N: {meta.get('adubacao','?')}\n")
        self.d_meta_box.insert("end", f"Nº de cortes (mow): {len(mow)}\n")
        if not df.empty:
            self.d_meta_box.insert(
                "end", f"Período observado: {df['date'].min().strftime('%d/%m/%Y')} "
                f"a {df['date'].max().strftime('%d/%m/%Y')}\n"
            )
        vars_disp = sorted(df["variable"].unique()) if not df.empty else []
        self.d_meta_box.insert("end", f"Variáveis disponíveis: {', '.join(vars_disp) if vars_disp else 'nenhuma'}")
        self.d_meta_box.configure(state="disabled")

        labels = [VAR_LABELS.get(v, v) for v in vars_disp]
        self.d_var_combo.configure(values=labels)
        if labels:
            preferred = "LAI (m²/m²)" if "LAI (m²/m²)" in labels else labels[0]
            self.d_var_combo.set(preferred)
        self._d_draw()

    def _d_draw(self):
        code = self._d_current_exp_code()
        info = self.dssat[code]
        df, mow = info["df"], info["mow"]
        var_label = self.d_var_combo.get()
        var_code = next((k for k, v in VAR_LABELS.items() if v == var_label), var_label)

        ax = self.d_ax
        ax.clear()
        ax.set_facecolor("#1f1f1f")
        if not df.empty and var_code in df["variable"].unique():
            sub = df[df.variable == var_code].sort_values("date")
            ax.plot(sub["date"], sub["value"], "-o", color="#5dade2", ms=4, lw=1.6)
        for md in mow:
            ax.axvline(md, color="#e74c3c", ls=":", lw=1, alpha=0.7)
        ax.set_title(f"{code} — {var_label}", color="white", fontsize=11)
        ax.set_ylabel(var_label, color="white")
        ax.tick_params(colors="white", labelrotation=25)
        for sp in ax.spines.values():
            sp.set_color("#555")
        self.d_fig.tight_layout()
        self.d_canvas.draw()

    # ============================================================ ABA 3 =====
    def _build_compare_tab(self, parent):
        ctk.CTkLabel(
            parent, text="Biomassa aérea acumulada (CWAD) — comparação entre todos os locais\n"
            "(CWAD é a única variável presente em TODOS os 5 experimentos)",
            font=ctk.CTkFont(size=13, weight="bold"),
        ).pack(pady=(14, 6))

        fig, axes = plt.subplots(len(self.dssat), 1, figsize=(9.5, 9), dpi=100, sharex=False)
        fig.patch.set_facecolor("#1a1a1a")
        if len(self.dssat) == 1:
            axes = [axes]
        colors = ["#5dade2", "#f39c12", "#2ecc71", "#e67e22", "#e74c3c"]
        for ax, (code, info), color in zip(axes, self.dssat.items(), colors):
            df, mow = info["df"], info["mow"]
            ax.set_facecolor("#1f1f1f")
            sub = df[df.variable == "CWAD"].sort_values("date") if not df.empty else pd.DataFrame()
            if not sub.empty:
                ax.plot(sub["date"], sub["value"], "-o", color=color, ms=3, lw=1.4)
            for md in mow:
                ax.axvline(md, color="#e74c3c", ls=":", lw=0.8, alpha=0.6)
            ax.set_ylabel("kg/ha", color="white", fontsize=8)
            ax.set_title(EXP_LABELS.get(code, code), color="white", fontsize=9, loc="left")
            ax.tick_params(colors="white", labelsize=7, labelrotation=20)
            for sp in ax.spines.values():
                sp.set_color("#555")
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)
        canvas.draw()

    # ============================================================ ABA 5 =====
    def _build_howmarkov_tab(self, parent):
        ctk.CTkLabel(
            parent, text="Por dentro do modelo de Markov: como a matriz é construída e como ela "
            "produz comportamento ao longo do tempo",
            font=ctk.CTkFont(size=14, weight="bold"),
        ).pack(pady=(10, 2))

        sub = ctk.CTkTabview(parent, width=1180, height=680)
        sub.pack(fill="both", expand=True, padx=8, pady=8)
        sA = sub.add("A · Construção da matriz")
        sB = sub.add("B · Diagrama de estados")
        sC = sub.add("C · Simulação (Monte Carlo)")

        self._build_construction_panel(sA)
        self._build_diagram_panel(sB)
        self._build_simulation_panel(sC)

    # ---- A: construção passo a passo da matriz a partir da sequência real -----
    def _build_construction_panel(self, parent):
        parent.grid_columnconfigure(0, weight=0)
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        side = ctk.CTkFrame(parent, width=300, corner_radius=0)
        side.grid(row=0, column=0, sticky="nswe")
        side.grid_propagate(False)

        ctk.CTkLabel(side, text="Genótipo", font=ctk.CTkFont(size=13, weight="bold")).pack(
            anchor="w", padx=14, pady=(14, 4)
        )
        self.c_geno = ctk.CTkSegmentedButton(
            side, values=["AGZG1219", "AGZG1501"], command=lambda v: self._c_reset()
        )
        self.c_geno.set("AGZG1219")
        self.c_geno.pack(fill="x", padx=14, pady=(0, 10))

        ctk.CTkLabel(
            side, text="A cada passo, pegamos um par de observações consecutivas "
            "(estado no dia t, estado no dia t+1) da sequência REAL de treino (2016) "
            "e somamos +1 na célula correspondente da tabela de contagens.",
            justify="left", wraplength=270, font=ctk.CTkFont(size=11), text_color="#b8c2cc",
        ).pack(fill="x", padx=14, pady=(0, 10))

        row_btns = ctk.CTkFrame(side, fg_color="transparent")
        row_btns.pack(fill="x", padx=14, pady=(0, 6))
        ctk.CTkButton(row_btns, text="🔄", width=40, command=self._c_reset).pack(side="left")
        self.c_play_btn = ctk.CTkButton(row_btns, text="▶ Reproduzir", command=self._c_toggle_play)
        self.c_play_btn.pack(side="left", padx=6)
        ctk.CTkButton(row_btns, text="Passo ➜", width=80, command=self._c_step).pack(side="left")

        self.c_pair_label = ctk.CTkLabel(side, text="", font=ctk.CTkFont(size=14, weight="bold"))
        self.c_pair_label.pack(anchor="w", padx=14, pady=(8, 2))
        self.c_progress_label = ctk.CTkLabel(side, text="", font=ctk.CTkFont(size=11))
        self.c_progress_label.pack(anchor="w", padx=14, pady=(0, 10))

        ctk.CTkLabel(side, text="Contagens acumuladas (linha = estado de origem)").pack(
            anchor="w", padx=14
        )
        self.c_counts_box = ctk.CTkTextbox(side, height=140, font=ctk.CTkFont(size=11, family="Courier"))
        self.c_counts_box.pack(fill="x", padx=14, pady=(2, 10))
        self.c_counts_box.configure(state="disabled")

        self.c_match_label = ctk.CTkLabel(side, text="", justify="left", wraplength=270,
                                           font=ctk.CTkFont(size=11))
        self.c_match_label.pack(fill="x", padx=14)

        main = ctk.CTkFrame(parent, fg_color="transparent")
        main.grid(row=0, column=1, sticky="nswe", padx=8, pady=8)
        main.grid_rowconfigure(0, weight=1)
        main.grid_columnconfigure(0, weight=1)
        self.c_fig, self.c_ax = plt.subplots(figsize=(7.2, 6.2), dpi=100)
        self.c_fig.patch.set_facecolor("#1a1a1a")
        self.c_canvas = FigureCanvasTkAgg(self.c_fig, master=main)
        self.c_canvas.get_tk_widget().grid(row=0, column=0, sticky="nswe")

        self.c_playing = False
        self._c_reset()

    def _c_reset(self):
        self.c_playing = False
        self.c_play_btn.configure(text="▶ Reproduzir")
        geno = self.c_geno.get()
        seq = (self.obs[(self.obs.experimento == geno) & (self.obs.amostra == "treino")]
               .sort_values("data")["estado"].tolist())
        self.c_seq = seq
        self.c_pairs = list(zip(seq[:-1], seq[1:]))
        self.c_pos = 0
        self.c_counts = np.zeros((5, 5))
        self._c_render()

    def _c_toggle_play(self):
        self.c_playing = not self.c_playing
        self.c_play_btn.configure(text="⏸ Pausar" if self.c_playing else "▶ Reproduzir")
        if self.c_playing:
            self._c_auto()

    def _c_auto(self):
        if not self.c_playing:
            return
        if self.c_pos >= len(self.c_pairs):
            self.c_playing = False
            self.c_play_btn.configure(text="▶ Reproduzir")
            return
        self._c_step()
        self.after(500, self._c_auto)

    def _c_step(self):
        if self.c_pos >= len(self.c_pairs):
            return
        a, b = self.c_pairs[self.c_pos]
        self.c_counts[a - 1, b - 1] += 1
        self.c_pos += 1
        self._c_render(last_pair=(a, b))

    def _c_render(self, last_pair=None):
        if last_pair:
            self.c_pair_label.configure(
                text=f"Transição contada:  Estado {last_pair[0]} → Estado {last_pair[1]}"
            )
        else:
            self.c_pair_label.configure(text="Clique em ▶ ou 'Passo' para começar")
        self.c_progress_label.configure(
            text=f"{self.c_pos} / {len(self.c_pairs)} transições da série de treino processadas"
        )

        self.c_counts_box.configure(state="normal")
        self.c_counts_box.delete("1.0", "end")
        header = "      " + "".join(f"E{j+1:<4}" for j in range(5))
        self.c_counts_box.insert("end", header + "\n")
        for i in range(5):
            linha = f"E{i+1}  | " + "".join(f"{int(self.c_counts[i,j]):<5}" for j in range(5))
            self.c_counts_box.insert("end", linha + "\n")
        self.c_counts_box.configure(state="disabled")

        # matriz normalizada parcial (o que já foi "aprendido" até este ponto)
        rows_sum = self.c_counts.sum(axis=1, keepdims=True)
        P_partial = np.divide(self.c_counts, rows_sum, out=np.zeros_like(self.c_counts),
                               where=rows_sum != 0)
        matrix_rows = [
            {"estado_origem": i + 1, "estado_destino": j + 1, "probabilidade": P_partial[i, j]}
            for i in range(5) for j in range(5) if P_partial[i, j] > 0
        ]
        highlight = last_pair[1] if last_pair else None
        draw_markov_diagram(
            self.c_ax, matrix_rows, highlight=highlight,
            title=f"Matriz parcial (aprendida até aqui) — {self.c_geno.get()}",
        )
        self.c_fig.tight_layout()
        self.c_canvas.draw()

        if self.c_pos >= len(self.c_pairs) and len(self.c_pairs) > 0:
            official = self.trans[self.trans.experimento == self.c_geno.get()]
            P_off = matrix_df_to_numpy(official)
            diff = float(np.abs(P_partial - P_off).max())
            self.c_match_label.configure(
                text=f"✅ Fim da série de treino. Diferença máxima vs. matriz oficial do "
                f"artigo: {diff:.3f} (esperado ≈ 0, pequenas diferenças vêm de arredondamento).",
                text_color="#7fd88f",
            )
        else:
            self.c_match_label.configure(text="")

    # ---- B: diagrama de estados completo, navegável ---------------------------
    def _build_diagram_panel(self, parent):
        parent.grid_columnconfigure(0, weight=0)
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        side = ctk.CTkFrame(parent, width=300, corner_radius=0)
        side.grid(row=0, column=0, sticky="nswe")
        side.grid_propagate(False)

        ctk.CTkLabel(side, text="Genótipo", font=ctk.CTkFont(size=13, weight="bold")).pack(
            anchor="w", padx=14, pady=(14, 4)
        )
        self.b_geno = ctk.CTkSegmentedButton(
            side, values=["AGZG1219", "AGZG1501"], command=lambda v: self._b_render()
        )
        self.b_geno.set("AGZG1219")
        self.b_geno.pack(fill="x", padx=14, pady=(0, 10))

        ctk.CTkLabel(side, text="Destacar estado de origem").pack(anchor="w", padx=14)
        self.b_state = ctk.CTkSegmentedButton(
            side, values=["Nenhum", "1", "2", "3", "4", "5"], command=lambda v: self._b_render()
        )
        self.b_state.set("Nenhum")
        self.b_state.pack(fill="x", padx=14, pady=(4, 14))

        ctk.CTkLabel(
            side, text="Cada círculo é um estado fenotípico. As setas mostram as "
            "probabilidades de transição estimadas na fase de treino (2016).\n\n"
            "Setas mais grossas = transições mais prováveis.\n"
            "Setas curvas voltando ao mesmo círculo = chance de permanecer no estado.\n\n"
            "Transições abaixo de 3% são omitidas para não poluir o desenho.",
            justify="left", wraplength=270, font=ctk.CTkFont(size=11), text_color="#b8c2cc",
        ).pack(fill="x", padx=14)

        main = ctk.CTkFrame(parent, fg_color="transparent")
        main.grid(row=0, column=1, sticky="nswe", padx=8, pady=8)
        main.grid_rowconfigure(0, weight=1)
        main.grid_columnconfigure(0, weight=1)
        self.b_fig, self.b_ax = plt.subplots(figsize=(7.2, 6.2), dpi=100)
        self.b_fig.patch.set_facecolor("#1a1a1a")
        self.b_canvas = FigureCanvasTkAgg(self.b_fig, master=main)
        self.b_canvas.get_tk_widget().grid(row=0, column=0, sticky="nswe")

        self._b_render()

    def _b_render(self):
        geno = self.b_geno.get()
        sub = self.trans[self.trans.experimento == geno]
        rows = sub.to_dict("records")
        hl = self.b_state.get()
        highlight = int(hl) if hl != "Nenhum" else None
        draw_markov_diagram(self.b_ax, rows, highlight=highlight,
                             title=f"Diagrama de transições — {geno} (treino 2016)")
        self.b_fig.tight_layout()
        self.b_canvas.draw()

    # ---- C: simulação de Monte Carlo ------------------------------------------
    def _build_simulation_panel(self, parent):
        parent.grid_columnconfigure(0, weight=0)
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        side = ctk.CTkFrame(parent, width=300, corner_radius=0)
        side.grid(row=0, column=0, sticky="nswe")
        side.grid_propagate(False)

        ctk.CTkLabel(side, text="Genótipo", font=ctk.CTkFont(size=13, weight="bold")).pack(
            anchor="w", padx=14, pady=(14, 4)
        )
        self.s_geno = ctk.CTkSegmentedButton(
            side, values=["AGZG1219", "AGZG1501"], command=lambda v: self._s_reset()
        )
        self.s_geno.set("AGZG1219")
        self.s_geno.pack(fill="x", padx=14, pady=(0, 10))

        row_btns = ctk.CTkFrame(side, fg_color="transparent")
        row_btns.pack(fill="x", padx=14, pady=(0, 6))
        ctk.CTkButton(row_btns, text="🔄", width=40, command=self._s_reset).pack(side="left")
        self.s_play_btn = ctk.CTkButton(row_btns, text="▶ Passo a passo", command=self._s_toggle_play)
        self.s_play_btn.pack(side="left", padx=6)
        ctk.CTkButton(side, text="⏩ Simular 2000 passos (rápido)", command=self._s_fast_forward).pack(
            fill="x", padx=14, pady=(6, 10)
        )

        self.s_draw_label = ctk.CTkLabel(side, text="", justify="left", wraplength=270,
                                          font=ctk.CTkFont(size=12))
        self.s_draw_label.pack(fill="x", padx=14, pady=(4, 6))

        ctk.CTkLabel(side, text="Régua de sorteio (nº aleatório entre 0 e 1)").pack(
            anchor="w", padx=14
        )
        self.s_ruler_fig, self.s_ruler_ax = plt.subplots(figsize=(3.0, 0.6), dpi=100)
        self.s_ruler_fig.patch.set_facecolor("#1a1a1a")
        self.s_ruler_canvas = FigureCanvasTkAgg(self.s_ruler_fig, master=side)
        self.s_ruler_canvas.get_tk_widget().pack(fill="x", padx=14, pady=(2, 12))

        ctk.CTkLabel(
            side, text="A régua mostra as faixas de probabilidade acumulada do estado atual. "
            "Sorteamos um número; a faixa em que ele cai define o próximo estado — "
            "é assim que a matriz 'produz' a trajetória.",
            justify="left", wraplength=270, font=ctk.CTkFont(size=11), text_color="#b8c2cc",
        ).pack(fill="x", padx=14)

        main = ctk.CTkFrame(parent, fg_color="transparent")
        main.grid(row=0, column=1, sticky="nswe", padx=8, pady=8)
        main.grid_rowconfigure(0, weight=1)
        main.grid_rowconfigure(1, weight=1)
        main.grid_columnconfigure(0, weight=1)
        self.s_fig1, self.s_ax1 = plt.subplots(figsize=(7.2, 3.0), dpi=100)
        self.s_fig1.patch.set_facecolor("#1a1a1a")
        self.s_canvas1 = FigureCanvasTkAgg(self.s_fig1, master=main)
        self.s_canvas1.get_tk_widget().grid(row=0, column=0, sticky="nswe")

        self.s_fig2, self.s_ax2 = plt.subplots(figsize=(7.2, 2.6), dpi=100)
        self.s_fig2.patch.set_facecolor("#1a1a1a")
        self.s_canvas2 = FigureCanvasTkAgg(self.s_fig2, master=main)
        self.s_canvas2.get_tk_widget().grid(row=1, column=0, sticky="nswe", pady=(8, 0))

        self.s_playing = False
        self._s_reset()

    def _s_reset(self):
        self.s_playing = False
        self.s_play_btn.configure(text="▶ Passo a passo")
        geno = self.s_geno.get()
        sub = self.trans[self.trans.experimento == geno]
        self.s_P = matrix_df_to_numpy(sub)
        self.s_state = 1
        self.s_trajectory = [1]
        self.s_visits = np.zeros(5)
        self.s_visits[0] = 1
        self._s_render(draw_info=None)

    def _s_toggle_play(self):
        self.s_playing = not self.s_playing
        self.s_play_btn.configure(text="⏸ Pausar" if self.s_playing else "▶ Passo a passo")
        if self.s_playing:
            self._s_auto()

    def _s_auto(self):
        if not self.s_playing:
            return
        if len(self.s_trajectory) >= 120:
            self.s_playing = False
            self.s_play_btn.configure(text="▶ Passo a passo")
            return
        self._s_step()
        self.after(500, self._s_auto)

    def _s_step(self):
        probs = self.s_P[self.s_state - 1]
        r = random.random()
        cum = np.cumsum(probs)
        nxt = int(np.searchsorted(cum, r) + 1)
        nxt = min(nxt, 5)
        self.s_state = nxt
        self.s_trajectory.append(nxt)
        self.s_visits[nxt - 1] += 1
        self._s_render(draw_info=(probs, r, nxt))

    def _s_fast_forward(self):
        for _ in range(2000):
            probs = self.s_P[self.s_state - 1]
            r = random.random()
            cum = np.cumsum(probs)
            nxt = int(np.searchsorted(cum, r) + 1)
            nxt = min(nxt, 5)
            self.s_state = nxt
            self.s_trajectory.append(nxt)
            self.s_visits[nxt - 1] += 1
        self._s_render(draw_info=None)

    def _s_render(self, draw_info):
        if draw_info:
            probs, r, nxt = draw_info
            self.s_draw_label.configure(
                text=f"Estava no Estado {self.s_trajectory[-2]}. Sorteado r = {r:.3f} → "
                f"caiu na faixa do Estado {nxt}."
            )
            self.s_ruler_ax.clear()
            self.s_ruler_ax.set_facecolor("#1a1a1a")
            self.s_ruler_ax.set_xlim(0, 1)
            self.s_ruler_ax.set_ylim(0, 1)
            self.s_ruler_ax.axis("off")
            cum = 0.0
            for k in range(5):
                w = probs[k]
                if w > 0:
                    self.s_ruler_ax.axvspan(cum, cum + w, color=STATE_COLORS[k + 1], alpha=0.85)
                cum += w
            self.s_ruler_ax.axvline(r, color="white", lw=2)
            self.s_ruler_fig.tight_layout(pad=0.2)
            self.s_ruler_canvas.draw()
        else:
            self.s_draw_label.configure(
                text="Estado inicial: 1 (pós-corte). Clique em 'Passo a passo' ou "
                "'Simular 2000 passos' para começar."
            )

        ax = self.s_ax1
        ax.clear()
        ax.set_facecolor("#1f1f1f")
        traj = self.s_trajectory[-150:]
        ax.step(range(len(traj)), traj, where="mid", color="#5dade2", lw=1.6)
        ax.scatter(range(len(traj)), traj, c=[STATE_COLORS[s] for s in traj], s=22, zorder=5)
        ax.set_yticks([1, 2, 3, 4, 5])
        ax.set_ylim(0.5, 5.5)
        ax.set_title(f"Trajetória simulada (últimos {len(traj)} passos) — {self.s_geno.get()}",
                     color="white", fontsize=10)
        ax.set_xlabel("passo simulado", color="white", fontsize=9)
        ax.tick_params(colors="white")
        for sp in ax.spines.values():
            sp.set_color("#555")
        self.s_fig1.tight_layout()
        self.s_canvas1.draw()

        ax2 = self.s_ax2
        ax2.clear()
        ax2.set_facecolor("#1f1f1f")
        freq = self.s_visits / self.s_visits.sum()
        stat = stationary_distribution(self.s_P)
        x = np.arange(1, 6)
        ax2.bar(x - 0.17, freq, width=0.34, color="#5dade2", label=f"Frequência simulada (n={int(self.s_visits.sum())})")
        ax2.bar(x + 0.17, stat, width=0.34, color="#f39c12", label="Distribuição estacionária (teórica)")
        ax2.set_xticks(x)
        ax2.set_title("Quanto tempo, no longo prazo, a planta passa em cada estado?",
                      color="white", fontsize=10)
        ax2.tick_params(colors="white")
        for sp in ax2.spines.values():
            sp.set_color("#555")
        ax2.legend(fontsize=7, facecolor="#2a2a2a", labelcolor="white")
        self.s_fig2.tight_layout()
        self.s_canvas2.draw()

    # ============================================================ ABA 4 =====
    def _build_info_tab(self, parent):
        box = ctk.CTkTextbox(parent, font=ctk.CTkFont(size=13), wrap="word")
        box.pack(fill="both", expand=True, padx=16, pady=16)
        box.insert("end", INFO_TEXT)
        box.configure(state="disabled")


INFO_TEXT = """COMO LER ESTE APLICATIVO
==========================

1) BOTÂNICA POR TRÁS DOS 5 ESTADOS FENOTÍPICOS

A alfafa rebrota a partir da coroa (crown), estrutura na base do caule que guarda
gemas dormentes protegidas do corte:

  Estado 1 (pós-corte): quase toda folha ativa foi removida; a planta depende de
  RESERVAS de carboidrato/N na raiz e coroa para reemitir hastes.

  Estado 2 (rebrota inicial): gemas da coroa emitem hastes novas com folhas
  trifolioladas pequenas; a planta ainda está "pagando" a reconstrução foliar.

  Estado 3 (crescimento intermediário): já autossustentável, expansão vegetativa
  acelerada.

  Estado 4 (fechamento de dossel): folhas de hastes vizinhas se sobrepõem, início
  do autossombreamento.

  Estado 5 (cobertura plena): platô fisiológico — ganho marginal de biomassa por
  folha extra cai bastante; na prática, corta-se por volta daqui.

2) O QUE MUDOU AO ANALISAR A BASE COMPLETA DO PROFESSOR (Alfalfa.zip)

  - AGZG1219 e AGZG1501 são AMBOS irrigados e adubados de forma semelhante
    (63 kg N/ha) — a ideia de "um é sequeiro, outro é irrigado" não se sustenta
    nos dados; foi um erro que corrigi durante a conversa.

  - Os códigos de talhão nos arquivos (ex.: "209", "209-a-b") aparecem tanto em
    AGZG1219 quanto em AGZG1501, sugerindo que os dois experimentos vêm de uma
    rede de talhões parcialmente sobreposta na Espanha, e não de duas cultivares
    geneticamente independentes.

  - AGZG1502 usa EXATAMENTE as mesmas observações de campo que AGZG1501 — a
    diferença entre os dois é apenas a regra de corte simulada no DSSAT
    (colheita observada / fixa a cada 40 dias / baseada em graus-dia). Isso
    confirma que este pacote de dados foi originalmente montado para comparar
    REGRAS DE MANEJO de corte no modelo CROPGRO-Forage, não genótipos.

  - AGZG1229 e ONCF1401 não têm LAI medido, só biomassa (CWAD) e herbagem (HERB).
    ONCF1401 é de Ottawa, Canadá — clima temperado frio, bem diferente da Espanha,
    útil para ilustrar o efeito do clima na dinâmica de corte-rebrota (aba 3).

3) RECOMENDAÇÃO PRÁTICA PARA O IX CBRG

  Antes de submeter/apresentar, alinhe com o Prof. Leonardo qual narrativa vocês
  vão sustentar:

  (a) tratar AGZG1219/AGZG1501 como "condições de manejo/local distintas" e
      reposicionar o artigo em torno da sensibilidade do modelo de Markov a
      diferentes regimes de corte e ambientes (defensável com os dados como estão), ou

  (b) manter a narrativa de "diferenciação genotípica", mas então será preciso
      uma fonte adicional que confirme que AGZG1219 e AGZG1501 correspondem a
      cultivares de alfafa geneticamente distintas (o pacote DSSAT usado aqui,
      por si só, não sustenta essa afirmação — ambos usam o código de cultivar
      genérico AL0001 no arquivo).

Qualquer uma das duas é um trabalho válido para o Eixo 2 do congresso — o
importante é a afirmação estar alinhada com o que os dados de fato mostram.
"""


if __name__ == "__main__":
    app = ExploradorAlfafa()
    app.mainloop()
