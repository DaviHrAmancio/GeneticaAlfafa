# Modelagem markoviana do LAI em alfafa — IX CBRG

Código, dados brutos, resultados e painel interativo da análise do **Índice de Área Foliar (LAI)** de alfafa (*Medicago sativa* L.) com cadeias de Markov, a partir de 65 observações DSSAT (Vale do Ebro, Aragón, Espanha). Resumo submetido ao IX Congresso Brasileiro de Recursos Genéticos.

Autores: Davi Henrique de Souza Amâncio; Leonardo Talavera Campos (IFC Campus Camboriú).

> **Transparência:** cada número pode ser rastreado até o arquivo bruto. Leia `docs/TRANSPARENCIA_E_LIMITACOES.md` antes de interpretar diferenças entre experimentos.

## Resultados-chave
| Métrica | Resultado |
|---|---|
| Observações de LAI (idênticas ao `LAID` bruto) | 65 |
| P12 AGZG1219 (treino) | 0,600 (3/5) com o CSV; 0,333 (2/6) com os limiares do resumo |
| P55 AGZG1501 (treino) | 0,333 (2/6) |
| Frobenius entre experimentos (treino) | 1,0650 (CSV); 1,2705 (limiares) |
| Diferença entre experimentos vs acaso | p = 0,78 (CSV); p = 0,62 (limiares) |
| Estabilidade 2016→2017 | não rejeitada (p = 0,17–0,69), com poucas transições de teste |
| Valor preditivo em 2017 | log-score −1,12 a −1,35 vs −1,61 do modelo uniforme |
| Transições para o estado 1 que são cortes (AGZG1219) | 4 de 4 |

## Como reproduzir
```bash
python -m pip install -r requirements.txt
python pipeline_alfafa_lai.py        # CSVs, figuras, snapshot, fontes do NotebookLM
python scripts/verify_output.py      # confere contra os números do resumo
python scripts/build_html.py         # painel interativo em html/painel_alfafa.html
```

## Estrutura
```
data/raw/            .ALT, .MOW, .ALX (DSSAT) e o CSV preparado pelo ChatGPT (apenas conferência)
pipeline_alfafa_lai.py
outputs/csv/         resultados em CSV
figures/             figuras 01–07
html/                painel interativo (arquivo único)
notebooklm_sources/  fontes numeradas + 00_NOTEBOOKLM_MASTER.md
docs/                metodologia, transparência, contexto, manifesto de evidências
scripts/             verificação e geração do HTML
config_snapshot.json semente, versões e hashes dos dados brutos
```

## Fontes
Dados de campo: arquivos DSSAT fornecidos pelo orientador. Contexto de mercado: ver `docs/CONTEXTO_ALFAFA.md`.
Finalidade acadêmica e educacional.
