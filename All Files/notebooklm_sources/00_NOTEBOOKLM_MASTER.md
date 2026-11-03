# MASTER DE DADOS PARA NOTEBOOKLM — Modelagem markoviana do LAI em alfafa (IX CBRG)

Repositório: https://github.com/DaviHrAmancio/GeneticaAlfafa

## Regras de integridade
Classificar cada conteúdo como: DADO EMPÍRICO (LAI observado, cortes) · PARÂMETRO ESTIMADO (matrizes, IC) ·
TESTE ESTATÍSTICO (permutação, bootstrap) · CONTEXTO DE LITERATURA (fontes públicas citadas). Não inventar números.
Os experimentos AGZG1219 e AGZG1501 são códigos de experimento DSSAT; confirmar com o orientador se representam genótipos distintos.

## Dado de origem
65 observações de LAI extraídas dos arquivos .ALT (campo LAID); cortes dos arquivos .MOW; experimentos na Espanha (Aragón/Vale do Ebro).
Treino = 2016; teste = 2017. Verificação do CSV contra o bruto:

| n_obs_bruto | n_obs_csv | n_pareadas | so_no_bruto | so_no_csv | max_dif_abs_LAI | amostra_divergente | estado_csv_diverge_dos_limiares |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 65 | 65 | 65 | 0 | 0 | 0.0000 | 0 | 9 |

## Duas regras de discretização (importante)
- `csv`: coluna `estado` do CSV preparado pelo ChatGPT. Reproduz exatamente os números do resumo submetido.
- `resumo`: limiares do resumo (1,0 / 2,5 / 3,8 / 5,0) aplicados ao LAI bruto. 9 das 65 observações mudam de estado.

## Números do resumo e sensibilidade

| regra_estados | metrica | valor | contagem | ic95_inf | ic95_sup |
| --- | --- | --- | --- | --- | --- |
| resumo | P12 AGZG1219 (treino) | 0.3333 | 2/6 | 0.0433 | 0.7772 |
| resumo | P55 AGZG1501 (treino) | 0.3333 | 2/6 | 0.0433 | 0.7772 |
| resumo | Frobenius entre experimentos (treino) | 1.2705 |  |  |  |
| csv | P12 AGZG1219 (treino) | 0.6000 | 3/5 | 0.1466 | 0.9473 |
| csv | P55 AGZG1501 (treino) | 0.3333 | 2/6 | 0.0433 | 0.7772 |
| csv | Frobenius entre experimentos (treino) | 1.0650 |  |  |  |

## Teste de diferença entre experimentos (permutação por ciclo de rebrota, treino)

| regra_estados | frobenius_observado | p_permutacao_por_ciclo | n_ciclos_treino | n_ciclos_AGZG1219 | n_ciclos_AGZG1501 | n_permutacoes |
| --- | --- | --- | --- | --- | --- | --- |
| resumo | 1.2705 | 0.6155 | 11 | 5 | 6 | 100000 |
| csv | 1.0650 | 0.7783 | 11 | 5 | 6 | 100000 |

## Estabilidade temporal 2016→2017 e capacidade preditiva

| regra_estados | experimento | n_transicoes_treino | n_transicoes_teste | linhas_comparaveis | frobenius_linhas_comparaveis | p_bootstrap_parametrico | logscore_matriz_propria | logscore_matriz_outro_experimento | logscore_matriz_agrupada | logscore_uniforme |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| resumo | AGZG1219 | 20 | 6 | 5 | 1.6708 | 0.1659 | -1.2435 | -1.3984 | -1.3193 | -1.6094 |
| resumo | AGZG1501 | 28 | 7 | 4 | 1.3506 | 0.3829 | -1.2955 | -1.3348 | -1.2451 | -1.6094 |
| csv | AGZG1219 | 20 | 6 | 5 | 1.3452 | 0.6877 | -1.1149 | -1.1683 | -1.0454 | -1.6094 |
| csv | AGZG1501 | 28 | 7 | 4 | 1.4317 | 0.2922 | -1.3473 | -1.4559 | -1.3517 | -1.6094 |

## Decomposição por cortes (a dinâmica é dominada pelo manejo de corte)

| regra_estados | experimento | amostra | n_transicoes | n_com_corte | intervalo_medio_dias | intervalo_mediano_dias | intervalo_max_dias | n_para_estado1 | n_para_estado1_por_corte | n_saidas_estado5 | n_saidas_estado5_por_corte | n_permanencias_estado5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| resumo | AGZG1219 | treino | 20 | 4 | 9.20 | 7.50 | 28 | 5 | 4 | 4 | 3 | 1 |
| resumo | AGZG1219 | teste | 6 | 1 | 8.33 | 7.50 | 18 | 2 | 1 | 1 | 1 | 0 |
| resumo | AGZG1501 | treino | 28 | 5 | 8.11 | 7.00 | 16 | 10 | 5 | 6 | 2 | 2 |
| resumo | AGZG1501 | teste | 7 | 1 | 13.00 | 7.00 | 49 | 2 | 1 | 3 | 1 | 2 |
| csv | AGZG1219 | treino | 20 | 4 | 9.20 | 7.50 | 28 | 4 | 4 | 4 | 3 | 1 |
| csv | AGZG1219 | teste | 6 | 1 | 8.33 | 7.50 | 18 | 1 | 1 | 1 | 1 | 0 |
| csv | AGZG1501 | treino | 28 | 5 | 8.11 | 7.00 | 16 | 7 | 5 | 6 | 2 | 2 |
| csv | AGZG1501 | teste | 7 | 1 | 13.00 | 7.00 | 49 | 2 | 1 | 3 | 1 | 2 |

## Métricas de rebrota por ciclo (sem discretização)

| variavel | experimento | n | media | dp | mediana | p_valor |
| --- | --- | --- | --- | --- | --- | --- |
| taxa_rebrota_LAI_dia | AGZG1219 | 4 | 0.1088 | 0.0233 | 0.1033 |  |
| taxa_rebrota_LAI_dia | AGZG1501 | 4 | 0.1520 | 0.0530 | 0.1526 |  |
| taxa_rebrota_LAI_dia | comparacao_Mann_Whitney | 8 |  |  |  | 0.3429 |
| LAI_max | AGZG1219 | 4 | 4.9950 | 0.8863 | 5.0150 |  |
| LAI_max | AGZG1501 | 4 | 5.5925 | 0.2351 | 5.7050 |  |
| LAI_max | comparacao_Mann_Whitney | 8 |  |  |  | 0.3429 |
| dias_ate_LAI_2.5 | AGZG1219 | 4 | 16.3627 | 3.3409 | 16.8853 |  |
| dias_ate_LAI_2.5 | AGZG1501 | 4 | 16.4896 | 5.6203 | 14.8371 |  |
| dias_ate_LAI_2.5 | comparacao_Mann_Whitney | 8 |  |  |  | 0.8857 |
| dias_ate_LAI_3.8 | AGZG1219 | 4 | 25.5314 | 6.2826 | 24.1461 |  |
| dias_ate_LAI_3.8 | AGZG1501 | 4 | 23.0341 | 9.0916 | 20.0903 |  |
| dias_ate_LAI_3.8 | comparacao_Mann_Whitney | 8 |  |  |  | 0.4857 |
| dias_ate_LAI_5.0 | AGZG1219 | 3 | 37.8584 | 8.2026 | 40.0000 |  |
| dias_ate_LAI_5.0 | AGZG1501 | 4 | 30.3666 | 12.5449 | 28.2071 |  |
| dias_ate_LAI_5.0 | comparacao_Mann_Whitney | 7 |  |  |  | 0.6286 |
| CWAD_kg_ha_no_corte | AGZG1219 | 4 | 3314.7500 | 839.8116 | 3365.0000 |  |
| CWAD_kg_ha_no_corte | AGZG1501 | 4 | 4485.0000 | 228.8376 | 4430.0000 |  |
| CWAD_kg_ha_no_corte | comparacao_Mann_Whitney | 8 |  |  |  | 0.0571 |

## Matrizes de transição (treino)

### AGZG1219 — treino 2016 — regra: Coluna `estado` do CSV (ChatGPT) — reproduz o resumo

Probabilidades:

| de/para | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| E1 | 0.000 | 0.600 | 0.200 | 0.200 | 0.000 |
| E2 | 0.000 | 0.000 | 0.667 | 0.333 | 0.000 |
| E3 | 0.000 | 0.000 | 0.000 | 0.500 | 0.500 |
| E4 | 0.167 | 0.000 | 0.000 | 0.500 | 0.333 |
| E5 | 0.750 | 0.000 | 0.000 | 0.000 | 0.250 |

Contagens:

| de/para | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| E1 | 0 | 3 | 1 | 1 | 0 |
| E2 | 0 | 0 | 2 | 1 | 0 |
| E3 | 0 | 0 | 0 | 1 | 1 |
| E4 | 1 | 0 | 0 | 3 | 2 |
| E5 | 3 | 0 | 0 | 0 | 1 |

### AGZG1501 — treino 2016 — regra: Coluna `estado` do CSV (ChatGPT) — reproduz o resumo

Probabilidades:

| de/para | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| E1 | 0.167 | 0.500 | 0.333 | 0.000 | 0.000 |
| E2 | 0.200 | 0.000 | 0.600 | 0.200 | 0.000 |
| E3 | 0.143 | 0.143 | 0.286 | 0.143 | 0.286 |
| E4 | 0.500 | 0.000 | 0.000 | 0.000 | 0.500 |
| E5 | 0.333 | 0.000 | 0.000 | 0.333 | 0.333 |

Contagens:

| de/para | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| E1 | 1 | 3 | 2 | 0 | 0 |
| E2 | 1 | 0 | 3 | 1 | 0 |
| E3 | 1 | 1 | 2 | 1 | 2 |
| E4 | 2 | 0 | 0 | 0 | 2 |
| E5 | 2 | 0 | 0 | 2 | 2 |

### AGZG1219 — treino 2016 — regra: Limiares do resumo aplicados ao LAI bruto

Probabilidades:

| de/para | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| E1 | 0.167 | 0.333 | 0.167 | 0.333 | 0.000 |
| E2 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 |
| E3 | 0.000 | 0.000 | 0.000 | 0.500 | 0.500 |
| E4 | 0.167 | 0.000 | 0.000 | 0.500 | 0.333 |
| E5 | 0.750 | 0.000 | 0.000 | 0.000 | 0.250 |

Contagens:

| de/para | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| E1 | 1 | 2 | 1 | 2 | 0 |
| E2 | 0 | 0 | 2 | 0 | 0 |
| E3 | 0 | 0 | 0 | 1 | 1 |
| E4 | 1 | 0 | 0 | 3 | 2 |
| E5 | 3 | 0 | 0 | 0 | 1 |

### AGZG1501 — treino 2016 — regra: Limiares do resumo aplicados ao LAI bruto

Probabilidades:

| de/para | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| E1 | 0.444 | 0.444 | 0.111 | 0.000 | 0.000 |
| E2 | 0.200 | 0.000 | 0.400 | 0.200 | 0.200 |
| E3 | 0.250 | 0.000 | 0.250 | 0.250 | 0.250 |
| E4 | 0.500 | 0.000 | 0.000 | 0.000 | 0.500 |
| E5 | 0.333 | 0.000 | 0.000 | 0.333 | 0.333 |

Contagens:

| de/para | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| E1 | 4 | 4 | 1 | 0 | 0 |
| E2 | 1 | 0 | 2 | 1 | 1 |
| E3 | 1 | 0 | 1 | 1 | 1 |
| E4 | 2 | 0 | 0 | 0 | 2 |
| E5 | 2 | 0 | 0 | 2 | 2 |


## Contexto da alfafa (fontes públicas)

| indicador | valor | unidade | ano | fonte |
| --- | --- | --- | --- | --- |
| Área de alfafa e misturas com gramíneas colhida para feno — EUA | 14.68 | milhões de acres (≈5,94 Mha) | 2025.0000 | Hay & Forage Grower, com dados USDA/NASS |
| Produção de feno de alfafa e misturas — EUA | 50.21 | milhões de toneladas (USDA usa short tons) | 2025.0000 | Hay & Forage Grower, com dados USDA/NASS |
| Área de alfafa — Argentina | 4000000 | ha | 2020.0000 | ANBA (cita área da América Latina) |
| Área de alfafa — Chile | 170000 | ha | 2020.0000 | ANBA |
| Área de alfafa — Peru | 120000 | ha | 2020.0000 | ANBA |
| Área de alfafa — Uruguai | 70000 | ha | 2020.0000 | ANBA |
| Área de alfafa — Brasil | 40000 | ha | 2020.0000 | ANBA; Embrapa (Tipificação de feno de alfafa para o Brasil) |
| Área de alfafa — Aragón (Espanha) | 57000 | ha (mais de) | 2020.0000 | SomosAlfalfa (matéria de jul/2020) |
| Exportações espanholas de alfafa desidratada | 750000 | t (mais de) | 2024.0000 | Diario del Campo |
| Comércio mundial de feno de alfafa — valor | 2.7 | bilhões US$ | 2016.0000 | ANBA |
| Participação de 5 importadores no consumo mundial de feno de alfafa | 84 | % | 2020.0000 | ANBA |
| Área somada de EUA, UE, Argentina, Rússia, Canadá e Austrália | 32 | milhões de ha | 2020.0000 | ANBA |
| Preço internacional do feno de alfafa (qualidade 'Good' a 'Supreme') | 220–280 | US$/t | 2023.0000 | Embrapa (Tipificação de feno de alfafa para o Brasil) |
| Potencial produtivo com manejo e irrigação | 25 | t MS/ha/ano |  | Embrapa (Tipificação de feno de alfafa para o Brasil) |
| Proteína bruta da alfafa | 15–22 | % |  | NAAIC (Importance of Alfalfa) |
| Modelo CROPGRO-Alfafa no Vale do Ebro — colheitas avaliadas | 58 | colheitas (990–4.617 kg/ha) | 2018.0000 | Malik et al., Agronomy Journal 110(5):1777-1790 |
