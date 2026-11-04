# Transparência e limitações

Este documento separa o que os dados sustentam do que o resumo submetido afirma. O resumo já foi submetido e não muda; a apresentação e o repositório podem ser mais precisos.

## 1. Natureza dos resultados
- **Empírico:** LAI observado (65 pontos), datas de corte, biomassa no corte.
- **Estimado:** matrizes de transição, ICs, distâncias de Frobenius.
- **Teste estatístico:** permutação por ciclo, bootstrap paramétrico.
- **Literatura:** números de contexto (área, comércio), todos com fonte e URL em `contexto_alfafa_fontes.csv`.

## 2. O que foi reproduzido
Com a coluna `estado` do CSV: **P12 (AGZG1219) = 0,600 (3/5)**, **P55 (AGZG1501) = 0,333 (2/6)** e **Frobenius = 1,0650**.

## 3. Limitações que afetam a interpretação

**3.1 Os "genótipos" são códigos de experimento.** AGZG1219 e AGZG1501 são nomes de experimentos DSSAT. Nos arquivos `.ALX` ambos usam o mesmo código de cultivar (`AL0001`), com nomes diferentes ("alfalfa" e "Aragon"), em manejos distintos: o 1219 é uma fazenda ("finca") com tratamentos de corte; o 1501 é uma parcela irrigada do centro de pesquisa de Aragón com 63 kg N/ha. Diferenças entre os dois misturam genótipo, local, manejo e irrigação. **Confirmar com o orientador se são genótipos distintos.**

**3.2 Amostra muito pequena.** Treino: 20 transições (1219) e 28 (1501). Teste: 6 e 7. P12 = 3/5 e P55 = 2/6 têm IC 95% de [0,15; 0,95] e [0,04; 0,78].

**3.3 A diferença entre experimentos não é distinguível do acaso.** Teste de permutação por ciclo (11 ciclos de treino): **p = 0,78** (regra `csv`) e **p = 0,62** (regra `resumo`). A afirmação "o modelo discrimina de forma precisa os padrões fisiológicos das cultivares" não é sustentada por este teste.

**3.4 P55 não difere entre os experimentos.** AGZG1219 tem P55 = 1/4 = 0,25 e AGZG1501 tem 2/6 = 0,33, com ICs quase idênticos. A "maior retenção de dossel" do 1501 tem diferença de uma observação.

**3.5 A regra de estados importa.** Com os limiares do resumo aplicados ao LAI bruto, **P12 (AGZG1219) cai de 0,600 para 0,333 (2/6)** e o Frobenius sobe para 1,2705. A coluna `estado` do CSV tem 9 de 65 estados que não seguem os limiares.

**3.5.1 Validação 2016→2017.** Os testes não rejeitam estabilidade (p entre 0,17 e 0,69), mas isso tem poder baixo com 6–7 transições. A distância treino-teste (1,35–1,67) é maior que a distância entre experimentos no treino. Por isso, "confirmou a estabilidade" deve ser apresentado como "não rejeitou".

**3.6 A dinâmica é dominada pelo manejo de corte.** No AGZG1219, **4 de 4** transições para o estado 1 coincidem com um corte (regra `csv`), e 3 de 4 saídas do estado 5 são cortes. P55 mede, em boa parte, a probabilidade de o próximo corte ainda não ter ocorrido. Excluindo as transições com corte, P55 do AGZG1501 sobe para 2/4. Os estados descrevem em grande parte o calendário de corte, não só a fisiologia da cultivar.

**3.7 Passos irregulares.** O intervalo entre observações varia de 1 a 49 dias. A cadeia trata cada par consecutivo como um passo, então "passos até o estado 5" não é tempo em dias.

**3.8 Dormência de inverno.** Os trechos de nov–fev têm LAI baixo e não são rebrotas comparáveis. Ficam fora das métricas por ciclo.

**3.9 Local único.** Todos os dados vêm de Aragón (Espanha). Os limiares de estados e as taxas de rebrota não são transferíveis ao Brasil sem calibração.

## 4. O que o repositório sustenta
- Os números do resumo são reproduzíveis a partir dos dados brutos.
- O LAI segue um ciclo em dente de serra ditado pelos cortes (figuras 01–02).
- A cadeia de Markov tem valor preditivo acima do modelo uniforme em 2017: log-score da matriz própria de −1,12 a −1,35 contra −1,61 do modelo uniforme, nos quatro casos (`estabilidade_temporal_2017.csv`).
- A matriz própria **não** supera de forma consistente a matriz do outro experimento nem a agrupada. Não há ganho de usar matrizes por experimento nesta amostra.
- Nos ciclos da estação principal (n = 4 por experimento), taxa de rebrota (0,109 vs 0,152 LAI/dia; p = 0,34) e dias até LAI 2,5 (16,4 vs 16,5) não diferem. A biomassa por corte é maior no 1501 (4.485 vs 3.315 kg/ha; p = 0,057, no limite), o que é compatível com a irrigação e a adubação desse experimento.

## 5. Sugestão de enquadramento para a apresentação
"Propomos um arcabouço reprodutível, de baixo custo, para transformar séries de LAI em indicadores probabilísticos. Com 65 observações, ele é demonstrativo: recupera a estrutura de dente de serra do manejo e tem valor preditivo acima do acaso, mas ainda não discrimina cultivares. O próximo passo é aumentar a amostra e controlar o manejo de corte."
