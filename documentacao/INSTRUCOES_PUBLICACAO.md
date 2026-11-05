# Como atualizar o repositório GeneticaAlfafa

Destino: https://github.com/DaviHrAmancio/GeneticaAlfafa

## Opção A — interface do GitHub (sem terminal)

1. Faça uma cópia de segurança do repositório atual (GitHub: **Code → Download ZIP**).
2. Extraia **este** arquivo ZIP; seus arquivos começam diretamente por `index.html`, `docs/` e `Apresentação/` (não existe uma pasta-contêiner extra).
3. No repositório `GeneticaAlfafa`, branch `main`, adicione/atualize os arquivos e pastas com o conteúdo extraído. O GitHub pode não substituir arquivos com mesmo nome apenas arrastando-os no uploader: se isso ocorrer, abra o arquivo atual, escolha **Edit** e substitua o conteúdo, ou use Git local (Opção B).
4. Faça um commit com a mensagem `Atualiza apresentação CBRG: Espanha + China e ajustes visuais`.
5. Em **Settings → Pages**, selecione `Deploy from a branch`, branch `main` e a pasta `/(root)` (o ZIP também contém um `docs/index.html` para quem usa `/docs`).
6. Após a publicação, confira `https://davihramancio.github.io/GeneticaAlfafa/` e teste os 17 slides, links e controles.

**Preserve** todos os arquivos antigos que não estão dentro deste ZIP; este pacote não é uma clonagem integral do repositório.

## Opção B — Git local (mais segura para atualizações)

```bash
git clone https://github.com/DaviHrAmancio/GeneticaAlfafa.git
cd GeneticaAlfafa
# Extraia/copiei aqui o conteúdo do ZIP, permitindo sobrescrever os arquivos coincidentes.
git status
git add index.html docs Apresentação dados documentacao .nojekyll README.md
git commit -m "Atualiza apresentação CBRG: Espanha + China e ajustes visuais"
git push origin main
```

A linha `git add` adiciona apenas os novos arquivos e pastas, sem apagar arquivos existentes. Confira `git status` antes de confirmar. Se o repositório já tiver um README importante, **mescle** o texto deste README em vez de substituí-lo cegamente.

## Publicação

- Nenhum `npm install`, servidor ou pipeline de build é necessário.
- A apresentação completa está no `index.html`; também é reproduzida em `docs/index.html` e `Apresentação/index.html` para compatibilidade com caminhos diferentes.
- O HTML é autocontido. Ao abrir em modo local, os gráficos funcionam sem conexão. Os links externos exigem internet.
- Este pacote não altera arquivos `.ALT` ou resultados originais do estudo espanhol. São preservados no repositório atual durante a mesclagem.
