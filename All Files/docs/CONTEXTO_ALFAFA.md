# Contexto: por que a alfafa importa

Todos os números abaixo vêm de fontes públicas, listadas com URL em `outputs/csv/contexto_alfafa_fontes.csv`. Anos e metodologias diferem entre fontes; a figura `07_contexto_area_alfafa.png` compara ordens de grandeza.

## Importância agronômica e econômica
- **Forrageira de alto valor nutritivo:** 15–22% de proteína bruta; uso principal em alimentação de vacas leiteiras (NAAIC).
- **Estados Unidos (2025):** 14,68 milhões de acres (≈5,94 Mha) de alfafa e misturas, 50,21 milhões de toneladas de feno e 3,42 t/acre (Hay & Forage Grower, com dados USDA/NASS).
- **Comércio internacional:** cerca de 6 milhões de t de feno de alfafa e US$ 2,7 bilhões em 2016. Cinco importadores (Arábia Saudita, Emirados, China, Japão, Coreia do Sul) concentram 84% do consumo (ANBA, 2020).
- **Preço internacional do feno (jul/2023):** US$ 220–280/t, conforme a qualidade (Embrapa).
- **Potencial produtivo:** até ~25 t de matéria seca/ha/ano com manejo e irrigação (Embrapa).

## Onde os dados foram coletados: Aragón, Espanha
- Os arquivos `.ALX` identificam os experimentos como "FINCA-SPAIN" (AGZG1219) e "Aragon, Montañana, Zaragoza, SPAIN" (AGZG1501). O 1501 cita o Centro de Investigación y Tecnología Agroalimentaria de Aragón (CITA).
- Esse conjunto vem do trabalho de adaptação do modelo CROPGRO à alfafa: Malik et al. (2018, *Agronomy Journal* 110:1777–1790) usaram parcelas irrigadas por aspersão no **Vale do Ebro**, com 58 colheitas (990–4.617 kg/ha) e LAI observado entre 0,1 e 6,7. O LAI do nosso CSV (0,17–6,51) está nessa faixa.
- **Aragón:** mais de 57.000 ha de alfafa, principalmente no Vale do Ebro, com cerca de 750.000 t (matéria de 2020). Cataluña tem 23.205 ha e ~250.000 t.
- **Espanha como exportadora:** as exportações de alfafa desidratada superaram 750.000 t em 2024 (alta de ~15%), com Arábia Saudita, Emirados e China entre os principais destinos; cerca de 2/3 em fardos e o resto em pellets (Diario del Campo). A mesma matéria de 2020 cita 74 fábricas de desidratação. *As fontes divergem sobre se a Espanha é a primeira ou a segunda exportadora mundial; por isso o número de ranking não foi usado.*
- **Relevância para o estudo:** na Espanha a alfafa é cultivo irrigado, de várias colheitas por ano, voltado à exportação. O manejo de corte (que domina a dinâmica do LAI neste estudo) é central para a produtividade e a qualidade do feno.

## Brasil
- Cerca de 40.000 ha de alfafa, principalmente no RS e PR; a Argentina tem ~4 milhões de ha (ANBA, 2020). Aragón (57.000 ha) tem uma área da mesma ordem de grandeza que o Brasil inteiro.
- O Brasil importa grandes volumes de sementes de alfafa da Califórnia e da Argentina (ANBA). A Embrapa aponta o potencial de expansão da cultura no Cerrado e na Caatinga.
- Este é um dos argumentos para estudar recursos genéticos e fenotipagem de alfafa: ampliar a base adaptada às condições brasileiras.

## Sobre dados de mercado financeiro (B3)
Não existe contrato futuro ou ativo de alfafa na B3, e não identifiquei empresa listada com exposição direta e relevante à alfafa. O ambiente de execução também não consegue acessar APIs de mercado (rede restrita). Em vez de forçar um proxy financeiro frágil, a introdução usa dados de área, produção e comércio de fontes públicas. Se quiser um proxy financeiro (por exemplo, preços de leite, milho ou soja), ele pode ser adicionado a partir de um script que você execute localmente, com a ressalva de que seria correlação, não medida de alfafa.

## Fontes
- [Hay & Forage Grower — Slow and steady hay gains in 2025](https://hayandforage.com/article-5581-Slow-and-steady-hay-gains-in-2025.html)
- [ANBA — Alfalfa in demand on international markets](https://anba.com.br/en/alfalfa-in-demand-on-international-markets/)
- [Embrapa — Tipificação de feno de alfafa para o Brasil](https://www.infoteca.cnptia.embrapa.br/infoteca/bitstream/doc/1165879/1/Tipificacao-de-feno-de-alfafa-para-o-Brasil.pdf)
- [SomosAlfalfa — Espanha, primeiro fabricante e exportador mundial de alfafa desidratada](https://somosalfalfa.com.ar/negocios-noticias/espana-afianzada-como-primer-fabricante-y-exportador-mundial-de-alfalfa-deshidratada/)
- [Diario del Campo — Exportações espanholas de alfafa em 2024](https://diariodelcampo.com/las-exportaciones-espanolas-de-alfalfa-aumentaron-un-quince-por-ciento-en-2024-con-arabia-emiratos-y-china-como-grandes-compradores/)
- [NAAIC — Importance of Alfalfa](https://naaic.org/Alfalfa/Importance.html)
- [Malik et al. (2018), Agronomy Journal](https://doi.org/10.2134/agronj2017.12.0680)
