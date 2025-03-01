# Pipeline prognóstico de Mielofibrose por evidência genética somática
- A mielofibrose é uma neoplasia de medula óssea cuja taxa de sobrevida pós-diagnóstico pode variar de meses à décadas.
- Este projeto auxilia a estimar o prognóstico da doença e sobrevida, por meio da filtragem de variantes somáticas.
- Foram filtradas 30 amostras do projeto "LMA Brasil" (WP048, WP093, WP087, WP060, WP056, WP066, WP064, WP072, WP078, WP285, WP280, WP274, WP276, WP270, WP216, WP306, WP297, WP291, WP295, WP204, WP160, WP164, WP162, WP212, WP170, WP196, WP180, WP188, WP140, WP126).

- Os arquivos VCF do projeto foram convertidos previamente da versão do genoma hg19 para hg38 utilizando o programa gatk LiftoverVcf com as posições hg19ToHg38.over.chain da UCSC, seguido de anotação pelo Ensembl-VEP (VEP annot).

- Este repositório compreende todos os arquivos essenciais para a filtragem e análise das variantes relacionadas ao prognóstico da Mielofibrose. Aqui, encontram-se os dados anotados para análise, o script de filtragem e outros arquivos relevantes.

- Os genes considerados para a filtragem foram obtidos do artigo fundador do GIPSS – Sistema de Prognóstico por Pontuação Inspirado em dados Genéticos (https://doi.org/10.1038/s41375-018-0107-z). São eles:
  - CALR
  - ASXL1
  - SRSF2
  - U2AF1
- Os prognósticos possíveis e sobrevida são:
  - Risco baixo (sobrevida média 26,4 anos)
  - Risco intermediário 1 (sobrevida média 10,3 anos)
  - Risco intermediário 2 (sobrevida média 4,6 anos)
  - Risco alto (sobrevida média 2,6 anos)
## (RESULTADOS)
- Os resultados e discussão do projeto podem ser acessados no seguinte caminho: <https://sites.google.com/view/g1-t5vsomticas/estudo-de-caso?authuser=0>.
- Abaixo segue o script utilizado para filtragem dos VCFs do estudo.

## 1. Preparação ambiente de trabalho
```
# Clonar github do projeto lmabrasil-hg48.git
%%bash
rm -rf lmabrasil-hg38
git clone https://github.com/renatopuga/lmabrasil-hg38
```
```
# Clonar github contendo amostras do projeto lmabrasil convertidos para hg38
!git clone https://github.com/bioinfoamos01/projetolma.git
```
```
# Remover arquivo README.md da pasta lma
!rm /content/projetolma/README.md
```

```
# Copiar os 30 arquivos da pasta projetolma(contendo amostras pós lift-over hg38) para a pasta /vep_output
!cp /content/projetolma/* /content/lmabrasil-hg38/vep_output
```

```
# Criar uma lista.txt com 4 genes de impacto para prognóstico de mielofibrose
!echo -e "CALR\nASXL1\nSRSF2\nU2AF1\n" > /content/lmabrasil-hg38/hpo/mielofibrose.txt
```

## 2. Instalação das Ferramentas Necessárias
### a) Instalação do BCFtools com plugin split-vep
O plugin permite extrair os campos de anotações estruturadas como INFO/CSQ criadas por bcftools/csq ou VEP (em nosso caso VEP).
Mais informações: https://samtools.github.io/bcftools/howtos/plugin.split-vep.html

```
%%bash
git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/samtools/bcftools.git
cd bcftools
make
make install
```

### b) Instalação do uDocker

Udocker é uma ferramenta essencial para executar containers Docker de forma simplificada em sistemas sem privilégios de root.

Neste workflow, sempre que empregamos o comando Udocker, o fazemos com a opção docker --allow-root.

A execução com privilégios de root é utilizada apenas temporariamente em nosso fluxo de trabalho pois o mesmo foi projetado para ser utilizado em um ambiente Google Colab; no entanto, essa prática não é recomendada.

Mais informações: [https://indigo-dc.github.io/udocker/](https://indigo-dc.github.io/udocker/)

```
%%bash
# Fonte: https://gist.github.com/mwufi/6718b30761cd109f9aff04c5144eb885
pip install udocker
udocker --allow-root install
```

### c) Download da imagem do ensembl-vep
Ensembl VEP é um conjunto de ferramentas usado para prever os impactos de variantes. Neste fluxo de trabalho, empregaremos o comando de filtragem do VEP para selecionar as variantes de interesse.

Devido à instalação do Udocker, é viável baixar a imagem do VEP utilizando o comando udocker --allow-root pull.

Mais informações: https://grch37.ensembl.org/info/docs/tools/vep/index.html

```
!udocker --allow-root pull ensemblorg/ensembl-vep
```

## 3. Filtragem
```
# Declarar lista com nome das 30 amostras do projeto lma brasil
SAMPLES = ["WP048","WP093","WP087","WP060","WP056","WP066","WP064","WP072","WP078","WP285","WP280","WP274","WP276","WP270","WP216","WP306","WP297","WP291","WP295","WP204","WP160","WP164","WP162","WP212","WP170","WP196","WP180","WP188","WP140","WP126"]
```

## **4. Análise**
Mais informações:  [https://pandas.pydata.org/](https://pandas.pydata.org/)

### a) Converter os outputs da filtragem em .csv

```

import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)

#Conversão -> output salvo na pasta /content/lmabrasil-hg38/csv_filtrados
!mkdir /content/lmabrasil-hg38/csv_filtrados
for i in SAMPLES:
  df = pd.read_csv(f'/content/lmabrasil-hg38/vep_output/liftOver_{i}_hg19ToHg38.vep.filter.tsv',sep='\t',index_col=False)
  df.to_csv(f'/content/lmabrasil-hg38/csv_filtrados/{i}_filtrado.csv', index=False)
```
### b) Gerar uma tabela bruta com o resultado da filtragem

```
import glob
import pandas as pd

# Lista pasta com os 30 arquivos csv(s)
csv_files = glob.glob('/content/lmabrasil-hg38/csv_filtrados/*.{}'.format('csv'))

#Une o resultado dos 30 arquivos csv(s) numa única tabela
df_concat = pd.concat([pd.read_csv(i) for i in csv_files], ignore_index=True)

#OUTPUT TABELA FINAL EM CSV
!mkdir /content/lmabrasil-hg38/tabela_final
df_concat.to_csv('/content/lmabrasil-hg38/tabela_final/tabela_final.csv', index=False)

df_concat
```
### c) Gerar uma tabela com os scores de risco

```
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)

df = pd.read_csv("/content/lmabrasil-hg38/tabela_final/tabela_final.csv")
#Extrai da tabela_final, em um csv resumido contendo informações das colunas "TumorID","HGVSc","SYMBOL", "HGVSp", "CHROM", "EXON"
df1 = df[["TumorID","SYMBOL", "HGVSc","HGVSp", "CHROM", "EXON"]]
df1.to_csv('/content/tabela_resumo_.csv', index=False)

#Converte a tabela resumida em uma lista
lista = df1.values.tolist()
S1 = []
variantes = []
var = []

#SEPARANDO O NOME DAS AMOSTRAS NÃO ENCONTRADAS NO FILTRO DE VARIANTES, NA VARIÁVEL S1, QUE COMPUTARÁ TODOS COMO +1 DEVIDO A AUSÊNCIA DE CALR TIPO1 PROTETOR)
for x in range(len(lista)):
  var = lista[x][0]
  variantes.append(var)

for i in SAMPLES:
  if i not in variantes:
    S1.append(i)

#dicionário para receber variantes, todos começam com pontuação +1, só são despontuados em -1 se for encontrado mutação tipo 1 em CALR
#PONTUA TODAS AS AMOSTRAS DA LISTA DE VARIANTES CONFORME CRITÉRIOS DA GISSP + GENES IDH1, IDH2 E EZH2
dicio = {}
for i in range(len(lista)):
  dicio[lista[i][0]] = 1


for i in range(len(lista)):
  if ("CALR" and "NM_004343.4:c.1099_1150del") in lista[i]: #mutação do tipo 1 em CARL
    for k, v in dicio.items():
      if lista[i][0] == k:
        dicio[lista[i][0]] +=-1

  if ("SRSF2" and "1/3") in lista[i]: #hostspot em exon 1 de SRSF2
    for k, v in dicio.items():
      if lista[i][0] == k:
        dicio[lista[i][0]] +=1

  if ("ASXL1" and "12/13") in lista[i]: #hotspot exon 12
    for k, v in dicio.items():
      if lista[i][0] == k:
        dicio[lista[i][0]] +=1

  if ("U2AF1" and "NM_006758.3:c.470A>C") in lista[i]: #hotspot mutação Q157 EM U2AF1
    for k, v in dicio.items():
      if lista[i][0] == k:
        dicio[lista[i][0]] +=1

for i in range(len(S1)):
  dicio[S1[i]] = 1


#Computando as pontuações usando um dicionário para armazenar as contas

contagem = {}

for i in dicio.values():
  if i in contagem:
    contagem[i] +=1
  else:
    contagem[i] =1

#verificação se a contagem possui resultados entre 0 à 3 ( 0 baixo, 1 inter-1, 2 inter2, 3-alto)
for i in range(4):
  if i not in contagem:
    contagem[i] = 0

```

### d) Gerar gráficos e tabelas com os resultados obtidos
```
import numpy as npimport
import matplotlib.pyplot as plt

#Prognóstico de Mielofibrose - variantes em genes de impacto no score
labels = 'CALR', 'ASLX1', 'U2AF1', 'SRSF2'
sections = [8, 5, 2, 2]

plt.pie(sections, labels=labels, autopct = '%1.1f%%')

plt.title('Frequência de variantes em genes de impacto para o prognóstico de mielofibrose (n=14)')
plt.show()
```
Na análise gráfica, dos 30 pacientes estudados, 14 apresentaram variantes nos genes-alvo selecionados para a avaliação prognóstica, conforme os critérios estabelecidos pelo GIPSS. Os 16 pacientes restantes, que não exibiram variantes, foram automaticamente atribuídos a uma pontuação de +1, uma vez que não possuíam a mutação protetora tipo 1 em CALR.
```
#Exibindo os dados usando um dataframe pandas

import pandas as pd

data = {"TOTAL DE PACIENTES": [contagem[0],contagem[1],contagem[2],contagem[3]],
        "SCORE": ["0", "1", "2", ">=3"],
        "RISCO": ["Baixo", "Intermediário 1", "Intermediário 2", "Alto"],
        "SOBREVIDA EM 5 ANOS (%)": ["94,0%", "73,0%", "40,0%", "14,0%"],
        "SOBREVIDA MÉDIA (ANOS)": ["26,4 anos", "10,3 anos", "4,6 anos", "2,6 anos"]}

dataf = pd.DataFrame(data)

dataf.style.set_caption('ESTRATIFICADOR DE PROGNÓSTICO DE MIELOFIBROSE BASEADADO EM GIPSS (sem cariotipagem)')
```

```

criterios = {"Penalidade": ["+1", "+1", "+1", "+1"],
        "Achados": ["Ausência de deleção de 52pb em CALR", "Mutação em Exon 1 de SRSF2", "Mutação em Exon 12 de ASXL1", "Mutação Q157 em U2AF1"]}
criteriosf = pd.DataFrame(criterios)
criteriosf.style.set_caption('CRITÉRIOS DE PONTUAÇÃO BASEADO EM GIPSS (adaptado)')
```

```
import numpy as npimport
import matplotlib.pyplot as plt

#Prognóstico de Mielofibrose - GIPSS adaptado
labels = 'BAIXO', 'INTERMEDIÁRIO 1', 'INTERMEDIÁRIO-2', 'ALTO'
sections = [6,21, 3, 0 ]

plt.pie(sections, labels=labels, autopct = '%1.1f%%')

plt.title('Frequência de prognóstico GIPSS* para Mielofibrose (n=30)')
plt.show()
```
Dos 30 pacientes analisados, todos receberam um escore de risco, variando de 0 a >= 3, após a filtragem de genes e variantes com impacto no prognóstico de mielofibrose (MF). Utilizamos o guideline do GIPSS para esta avaliação, porém sem considerar a avaliação do cariótipo.

### Autores 📃
Amós Eduardo - codificação de script\
Renato puga - codificação de script\
Rabiana Rocha - revisão científica\
Luiz Gustavo - revisão científica\
Thiago adalton - front-end (google sites)\
Clara Sabença - front-end (google sites)

