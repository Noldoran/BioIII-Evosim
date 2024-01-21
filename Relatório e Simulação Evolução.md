# Introdução

Ao atuar sobre mutações condicionalmente benéficas, processos de seleção podem criar frequências alélicas pouco esperadas dentro de uma população. Um tal caso não-trivial é aquele na qual a viabilidade de um dado genótipo varia conforme o tempo.
Nesse trabalho, buscamos implementar tal caso em uma simulação simples que modela a evolução de uma população ao longo de um número de gerações.
Para esse fim, os seguintes parâmetros serão utilizados como entrada para a simulação:
- Frequência inicial $f_{\text{A}}$ do alelo $\text{A}$ 
- Tamanho da população $N$
- Funções da viabilidade de cada genótipo dado um número $g$ de gerações decorrido:
	- $w_{\text{AA}}(g)$ função que descreve a viabilidade de $\text{AA}$ na geração $g$
	- $w_{\text{Aa}}(g)$ função que descreve a viabilidade de $\text{Aa}$ na geração $g$
	- $w_{\text{aa}}(g)$ função que descreve a viabilidade de $\text{aa}$ na geração $g$

# Modelos determinísticos

Para que possamos validar a simulação, podemos comparar seus resultados com um modelo determinístico. Para esse fim, usaremos uma modificação do equilíbrio de Hardy-Weinberg que leva em conta a adaptabilidade relativa dos três possíveis genótipos $\text{AA}$, $\text{Aa}$ e $\text{aa}$.
$$\overline{w} = f_{\text A}^2 \cdot w_{\text{AA}} + 2f_{\text A}f_{\text a} \cdot w_{\text{Aa}} + f_{\text a}^2 \cdot w_{\text{aa}}$$
$$f'_{\text A} = \frac{f_{\text A}^2 \cdot w_{\text{AA}} + f_{\text A}f_{\text a} \cdot w_{\text{Aa}}}{\overline{w}}$$

Onde $\overline{w}$ é a fitness média da população e $f_{\text A}'$ é a frequência do alelo $\text{A}$ na geração seguinte.

Vamos também avaliar a evolução da distribuição genotípica da população. Dessa forma, temos os seguintes valores para a frequência de cada genótipo na geração seguinte:
$$\begin{align}
f_{\text{AA}}' &= f_{\text A}^2 \cdot w_{\text{AA}} \\
f_{\text{Aa}}' &= 2f_{\text A}f_{\text a} \cdot w_{\text{Aa}} \\
f_{\text{aa}}' &= f_{\text a}^2 \cdot w_{\text{aa}} \\
\end{align}$$

Crucialmente, deve-se observar que, com valores adaptativos variáveis, essas quantidades passam a ser as funções $w_{\text{AA}}(g)$, $w_{\text{Aa}}(g)$, e  $w_{\text{aa}}(g)$, dependentes de $g$, o número de gerações passadas - aqui usado como unidade discreta de tempo.

Dessa forma, as expressões usadas para a frequência do alelo $\text{A}$ e frequência de cada genótipo serão, respectivamente
$$\begin{align}
f'_{\text A} &= \frac{f_{\text A}^2 \cdot w_{\text{AA}}(g) + f_{\text A}f_{\text a} \cdot w_{\text{Aa}}(g)}{f_{\text A}^2 \cdot w_{\text{AA}}(g) + 2f_{\text A}f_{\text a} \cdot w_{\text{Aa}}(g) + f_{\text a}^2 \cdot w_{\text{aa}}(g)}
\end{align}$$
$$\begin{cases}
f_{\text{AA}}' &= f_{\text A}^2 \cdot w_{\text{AA}}(g) \\
f_{\text{Aa}}' &= 2f_{\text A}f_{\text a} \cdot w_{\text{Aa}}(g) \\
f_{\text{aa}}' &= f_{\text a}^2 \cdot w_{\text{aa}}(g) \\
\end{cases}$$

# Método de simulação

A simulação consiste de uma simples seleção aleatória de pares reprodutivos, enviesada pelo valor adaptativo de cada indivíduo. Para cada indivíduo a ser criado, selecionamos seus progenitores da seguinte forma na função `simulate_generation`:

Dada uma lista de genótipos `["AA", "Aa", "aa"]` calcula-se, por exemplo, a probabilidade de reprodução do genótipo $\text{AA}$ na geração `g` como `freq["AA"][g] * W["AA"](g)`. Uma vez que o par foi selecionado, chama-se a função `simulate_sex`, que escolhe aleatoriamente e sem viés um alelo de cada progenitor que formarão o genótipo da prole. Essa prole então é adicionada a população da nova geração `g+1`.

Esse processo é repetido até que a nova geração tenha $N$ integrantes.

O processo de simulação de uma geração, em turno, é repetido um número `G` de gerações, para que tenhamos os dados completos da simulação. Uma simulação de $10^4$ indivíduos ao longo de $100$ gerações (`G = 100`, `N = 10000`) leva cerca de sete segundos para completar. Há possibilidade de melhoria na eficiência utilizando-se a biblioteca `numpy`, porém nós não exploramos essa possibilidade em favor da clareza do código.

# Exemplo: Anemia Falciforme

Um caso bastante estudado de adaptabilidade variável de um traço *single-locus* é a interação entre a mutação do Traço Falciforme e a pressão ambiental causada pela Malária.

Na falta da presença da malária no ambiente, o indivíduo sem a mutação ($\text{AA}$), geralmente tem *fitness* mais alto do que aqueles que tem uma cópia do alelo mutado ($\text{Aa}$). Ainda mais, o indivíduo com duas cópias desse alelo ($\text{aa}$) desenvolve uma condição debilitante chamada de Anemia Falciforme. Porém, na presença de malária, indivíduos heterozigotos demonstram uma taxa muito reduzida de mortalidade devida à malária. 

Vale notar que existe também uma segunda variação do alelo do traço falciforme que aparenta ser desfavorecido em relação à mutação estudada aqui, conforme discutido em [P W Hedrick (2011)](https://doi.org/10.1038%2Fhdy.2011.16). Neste trabalho usaremos como exemplo o alelo denominado $\text{C}$ no artigo supracitado.

Para começar, utilizaremos os valores adaptativos encontrados em [P W Hedrick (2011)](https://doi.org/10.1038%2Fhdy.2011.16) para cada um dos três genótipos em Burquina Fasso comparando também com a frequência da mutação encontrada nessa região. Além disso, os seguintes valores relativos de fitness:

|  | $\text{AA}$ | $\text{Aa}$ | $\text{aa}$ |
| ---- | ---- | ---- | ---- |
| Malária ausente | 1 | 1 | 0.95 |
| Malária presente | 0.861 | 0.935 | 1 |

Os valores para o ambiente de *malária ausente* foram escolhidos arbitrariamente, e seu ajuste pode gerar resultados interessantes por exemplo no caso em que a pressão seletiva na ausência da malária seja forte o bastante para  eliminar completamente o alelo $\text{a}$ antes da introdução da malária na população.

Abaixo estão os resultados de uma simulação com os parâmetros de fitness listados a cima, bem como $N = 10^5$, $f_{\text{A}} = 0.75$ ao longo de $100$ gerações.
![Fig.1](https://raw.githubusercontent.com/Noldoran/BioIII-Evosim/main/allele_freq.png)
![Fig.2](https://raw.githubusercontent.com/Noldoran/BioIII-Evosim/main/genotype_distribution.png)

Note a linha vertical pontilhada, que denota a geração na qual a malária foi introduzida na população (na simulação mostrada, a geração $g = 50$)

Fica claro que a simulação, que ocorre de forma discretizada e estocástica, apresenta resultados muito próximos do modelo teórico determinístico apresentado.


# Resultados e conclusão

Pudemos desenvolver uma simulação que concorda com modelos determinísticos descritos na literatura. A simulação aparenta ser bastante robusta a uma ampla variação dos parâmetros de entrada, e não encontramos nenhuma tal combinação de parâmetros que cause uma deviação inesperada do modelo teórico.

Uma possibilidade interessante de se explorar e provavelmente mais realista seria o uso de funções de fitness $w(g)$ com característica contínua. A simulação comporta essa possibilidade, porem o exemplo apresentado faz uso de funções que não apresentam continuidade.

O código da simulação, bem como geração dos gráficos apresentados e este próprio relatório podem ser encontrados no repositório https://github.com/Noldoran/BioIII-Evosim.

Há ainda questões a se explorar que, apesar de caírem fora do escopo do presente trabalho, merecem ser exploradas com mais profundidade. Por exemplo, a concordância entre a simulação e a distribuição genotípica encontrada no mundo real, ou a expansão para acomodar mais de duas variações de cada alelo, como é o caso das variações $\text{HbS}$ e $\text{HbC}$ do traço falciforme.