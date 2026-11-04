# Manifesto de evidências

Mapeia cada afirmação do resumo e do repositório ao arquivo que a sustenta.

| Afirmação | Valor | Arquivo | Natureza |
|---|---|---|---|
| 65 observações de LAI, idênticas ao `LAID` bruto | 65/65, diferença máx. 0 | `verificacao_csv_vs_bruto.csv` | Empírico |
| 9 estados do CSV não seguem os limiares do resumo | 9 de 65 | `verificacao_csv_vs_bruto.csv`, `observacoes_lai.csv` | Empírico |
| P12 AGZG1219 (treino) | 0,600 (3/5) com CSV; 0,333 (2/6) com limiares | `sensibilidade_regra_de_estados.csv` | Estimado |
| P55 AGZG1501 (treino) | 0,333 (2/6) nas duas regras | `sensibilidade_regra_de_estados.csv` | Estimado |
| Frobenius entre experimentos (treino) | 1,0650 (CSV); 1,2705 (limiares) | `distancias_frobenius.csv` | Estimado |
| Diferença entre experimentos vs acaso | p = 0,78 (CSV); p = 0,62 (limiares) | `teste_permutacao_entre_experimentos.csv` | Teste |
| Estabilidade 2016→2017 | p = 0,17–0,69; 6–7 transições de teste | `estabilidade_temporal_2017.csv` | Teste |
| Valor preditivo acima do uniforme | −1,12 a −1,35 vs −1,61 | `estabilidade_temporal_2017.csv` | Estimado |
| Transições para o estado 1 são cortes | 4/4 (AGZG1219, CSV) | `decomposicao_cortes.csv` | Empírico |
| Rebrota por ciclo | 0,109 vs 0,152 LAI/dia (p = 0,34) | `resumo_ciclos.csv`, `ciclos_rebrota.csv` | Estimado |
| Biomassa por corte | 3.315 vs 4.485 kg/ha (p = 0,057) | `resumo_ciclos.csv` | Empírico |
| Distribuição estacionária e passos até o estado 5 | por experimento | `propriedades_cadeia.csv` | Estimado |
| Local dos experimentos | Aragón, Espanha | `data/raw/*.ALX` | Empírico |
| Números de contexto (área, comércio, preços) | com fonte e URL | `contexto_alfafa_fontes.csv` | Literatura |

## Figuras
| Figura | Mostra | Arquivos de origem |
|---|---|---|
| `01_serie_LAI_cortes.png` | LAI observado e cortes | `observacoes_lai.csv`, `data/raw/*.MOW` |
| `02_trajetorias_por_ciclo.png` | Rebrota por ciclo | `observacoes_lai.csv` |
| `03_matrizes_treino_2016.png` | Matrizes com contagens, duas regras | `matrizes_transicao_long.csv` |
| `04_incerteza_P12_P55.png` | ICs de P12 e P55 | `matrizes_transicao_long.csv` |
| `05_permutacao_frobenius.png` | Nulo da distância entre experimentos | `teste_permutacao_entre_experimentos.csv` |
| `06_validacao_preditiva_2017.png` | Log-score em 2017 | `estabilidade_temporal_2017.csv` |
| `07_contexto_area_alfafa.png` | Área de alfafa por país/região | `contexto_alfafa_fontes.csv` |

Regeneração: `python pipeline_alfafa_lai.py` e verificação: `python scripts/verify_output.py`.
