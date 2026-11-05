# Genética da Alfafa | Apresentação expandida para o CBRG

**Apresentação interativa (17 slides)**: modelagem estocástica markoviana do índice de área foliar (LAI), contextualização da importância agronômica da alfafa e análise descritiva das bases experimentais chinesas apresentadas em Yang et al. (2024).

## Abrir a apresentação

- **Pela página principal do repositório / GitHub Pages:** [index.html](index.html)
- **Pela pasta da apresentação:** [Apresentação/index.html](Apresenta%C3%A7%C3%A3o/index.html)
- **Com GitHub Pages configurado para `/docs`:** [docs/index.html](docs/index.html)
- Navegação: setas do teclado, barra inferior, botões e filtros interativos.
- Funciona offline, sem precisar de Node, servidor, build, bibliotecas ou fontes externas. Os links para artigos e repositórios abrem páginas externas e precisam de internet.

## Conteúdo do pacote

| Caminho | Função |
| --- | --- |
| `index.html` | Entrada principal para GitHub Pages (branch `main`, diretório `/(root)`) |
| `docs/index.html` | Entrada alternativa se GitHub Pages usa branch `main`, diretório `/docs` |
| `Apresentação/index.html` | Apresentação dentro da pasta usada pelo projeto |
| `Apresentação/Alfafa_ importância e modelagem markoviana do LAI.html` | Arquivo com o nome original, atualizado para preservar possíveis links antigos |
| `dados/China/leaf-area-and-leaf-area-index-main.zip` | Arquivo de dados original enviado para esta atualização (planilhas e subarquivos do estudo chinês) |
| `documentacao/Notas_Auditoria_Alfafa_CBRG.txt` | Notas de proveniência e ressalvas científicas |
| `documentacao/INSTRUCOES_PUBLICACAO.md` | Procedimento de atualização do GitHub |
| `.nojekyll` | Evita processamento Jekyll desnecessário no GitHub Pages |

## Bases científicas e escopo

1. **Espanha**: análise markoviana original embutida no HTML e dados já disponíveis no repositório de origem: [GeneticaAlfafa](https://github.com/DaviHrAmancio/GeneticaAlfafa).
2. **China**: material enviado pelo pesquisador, associado a [Yang et al. (2024), Frontiers in Plant Science](https://pmc.ncbi.nlm.nih.gov/articles/PMC10896985/) e ao [repositório dos dados](https://github.com/im1robot/leaf-area-and-leaf-area-index).
3. **Estatística**: as matrizes, validação e testes markovianos correspondem ao estudo espanhol. Os gráficos chineses são descritivos; não são testes novos nem uma validação markoviana já realizada. A estação climática de Zaragoza é referência regional, não o registro meteorológico exato da parcela espanhola.

## Atenção ao atualizar

Este ZIP é um **pacote completo da apresentação e das fontes chinesas fornecidas**, não é uma cópia integral de todos os arquivos preexistentes do repositório público (como `Resumo Submetido`, `old 20K` e os arquivos brutos `.ALT`). **Não exclua essas pastas nem esses dados**. Faça upload/mesclagem dos arquivos deste ZIP sobre o repositório existente, substituindo os de mesmo nome e preservando os outros.

Veja [instruções de publicação](documentacao/INSTRUCOES_PUBLICACAO.md).
