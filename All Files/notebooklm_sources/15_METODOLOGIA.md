# Metodologia

## Dados de origem
- **Observado (`.ALT`)**: campo `LAID` (LAI, m²/m²), além de `CWAD` (biomassa aérea, kg/ha) e `CHTD`. Formato DSSAT; `-99` = ausente.
- **Cortes (`.MOW`)**: datas dos eventos de corte simulados/observados no CROPGRO-Perennial Forage.
- **Experimento (`.ALX`)**: local, solo, manejo e cultivar.
- Datas DSSAT `AADDD` (ex.: `16118` = dia 118 de 2016 = 27/04/2016).
- O CSV `Alfafa_LAI_Markov_treino_teste.csv` foi preparado pelo ChatGPT. O pipeline confere que seus 65 valores de LAI são **idênticos** ao `LAID` bruto (diferença máxima = 0). A coluna `estado` é a única parte derivada.

## Divisão temporal
- **Treino:** 2016. **Teste fora da amostra:** 2017 (mesma divisão do resumo).
- Transições só são formadas entre observações consecutivas do **mesmo experimento e mesma amostra** (nenhuma transição atravessa 2016→2017).

## Estados do LAI (duas regras, ambas reportadas)
| Estado | Nome | Limiar do resumo |
|---|---|---|
| 1 | Pós-corte | LAI < 1,0 |
| 2 | Rebrota inicial | 1,0 ≤ LAI < 2,5 |
| 3 | Crescimento intermediário | 2,5 ≤ LAI < 3,8 |
| 4 | Fechamento de dossel | 3,8 ≤ LAI < 5,0 |
| 5 | Cobertura plena | LAI ≥ 5,0 |

- **Regra `csv`**: usa a coluna `estado` do CSV. É a que reproduz exatamente os números do resumo submetido.
- **Regra `resumo`**: aplica os limiares acima ao `LAID` bruto. Em 9 das 65 observações o estado difere do CSV (valores logo abaixo de um limiar, como LAI = 0,90 e 2,45).

## Cadeia de Markov
Tempo discreto, ordem 1, estimada por frequência relativa: `P_ij = n_ij / n_i·`. O "passo" é o intervalo entre observações consecutivas (1 a 49 dias; mediana ≈ 7 dias), **não** um passo de tempo fixo.
Linhas sem observação ficam com zero na matriz (convenção do resumo, necessária para reproduzir a distância de Frobenius).

## Distância de Frobenius
`F(A,B) = ||A − B||_F`. Reproduz 1,0650 (treino, regra `csv`).

## Incerteza e testes (adicionados ao resumo)
1. **IC 95% de Clopper-Pearson** para cada probabilidade de transição.
2. **Teste de permutação por ciclo de rebrota** para a distância entre experimentos (treino). Cada ciclo entre dois cortes é a unidade permutada, porque transições dentro de um ciclo são dependentes. 100.000 permutações.
3. **Estabilidade temporal 2016→2017**: bootstrap paramétrico, 100.000 réplicas (a matriz de 2016 gera os dados de 2017 com os mesmos tamanhos de linha) nas linhas comparáveis; mais um **escore preditivo** (log-verossimilhança média por transição, matrizes suavizadas com α = 0,5) contra o modelo uniforme, a matriz do outro experimento e a matriz agrupada.
4. **Propriedades da cadeia**: distribuição estacionária, permanência e número médio de passos até o estado 5 (somente quando a cadeia é irredutível).

## Análise por ciclo de corte (sem discretizar)
Cada ciclo é um segmento entre dois cortes do mesmo ano. Por ciclo: LAI mínimo, LAI máximo, taxa de rebrota `(máx − mín)/dias`, dias até LAI 2,5 / 3,8 / 5,0 (interpolação linear) e biomassa `CWAD` no corte final.
Entram nas comparações apenas ciclos da estação principal (rebrota iniciada de abril a outubro), encerrados por corte e com amplitude ≥ 1 de LAI. Isso exclui os trechos de dormência de inverno.
Observação no dia do corte: é tratada como pré-corte quando a observação seguinte tem LAI menor (caso do AGZG1501), e como pós-corte caso contrário (caso do AGZG1219).

## Reprodutibilidade
- `python pipeline_alfafa_lai.py` regenera todos os CSVs, figuras e o `config_snapshot.json`.
- Semente fixa: `20261103`. Os hashes SHA-256 dos arquivos brutos ficam no snapshot.
- `python scripts/verify_output.py` confere os resultados-chave.
