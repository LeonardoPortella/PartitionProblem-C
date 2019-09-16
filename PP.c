#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include </home/leonardo/IA/PP/Util.c>
#include </home/leonardo/IA/PP/GeraTeste.c>

const int tamPopulacao = 1000;
const int percSobreviventes = 10;
const int percMutacao = 2;
const int repeticoes = 1000;
const int percElite = 1;
const int percRenovacao = 20;
const int limSemMelhorarFit = 500;
	
//============================================================================================================

typedef struct 
{
	int *cromossomo;
	int aptidao;
	int crossover;
	int elite;
	float aptidaoRoleta;

} Cromossomo;

//============================================================================================================

int ordena(const void *a, const void *b) 
{
	return ( a - b );
}

//============================================================================================================

int ordenaAptidao(const void *a, const void *b) 
{
	Cromossomo *cromossomo0 = (Cromossomo *)a;
	Cromossomo *cromossomo1 = (Cromossomo *)b;

	return ( cromossomo0->aptidao - cromossomo1->aptidao );
}

//============================================================================================================

int ordenaCrossover(const void *a, const void *b) 
{
	Cromossomo *cromossomo0 = (Cromossomo *)a;
	Cromossomo *cromossomo1 = (Cromossomo *)b;

	return ( cromossomo1->crossover - cromossomo0->crossover ); //Desc
}

//============================================================================================================

void PrintSolucao(int *cromossomo, long int *genes, int qtdGenes)
{
	int *solucao0 = (int *)malloc(sizeof(int));
	int tamSol0 = 0;
	long int peso0 = 0;
	int *solucao1 = (int *)malloc(sizeof(int));;
	int tamSol1 = 0;
	long int peso1 = 0;
	long int pesoDif = 0;
	
	for(int i = 0; i < qtdGenes; ++i)
	{
		if(cromossomo[i] == 0)
		{		
			solucao0 = (int *)realloc(solucao0, ++tamSol0 * sizeof(int));
			solucao0[tamSol0 - 1] = genes[i];  
			peso0 += genes[i]; 
		}
		else
		{
			solucao1 = (int *)realloc(solucao1, ++tamSol1 * sizeof(int));
			solucao1[tamSol1 - 1] = genes[i];  
			peso1 += genes[i];  
		}
	}

	printf("Peso %li ", peso0);
	
	if(qtdGenes <= 15)
		PrintLista(solucao0, tamSol0);
	
	printf("Peso %li ", peso1);
	
	if(qtdGenes <= 15)
		PrintLista(solucao1, tamSol1);

	pesoDif = peso0 - peso1;

	if(pesoDif < 0)
		pesoDif *= -1;

	printf("Diferenca Peso %li\n", pesoDif);
	
	free(solucao0);
	free(solucao1);
}		

//============================================================================================================

void printPopulacao(Cromossomo *populacao, int tamPopulacao, int qtdGenes)
{
	for(int i = 0; i < tamPopulacao; ++i)
	{
		printf("cromossomo ");
		PrintLista(populacao[i].cromossomo, qtdGenes);
		printf("apt %i\n", populacao[i].aptidao);
	}
}

//============================================================================================================

void CromossomoAleatorio(Cromossomo *individuo, int qtdGenes)
{
	individuo->aptidao = -1;
	individuo->aptidaoRoleta = -1;
	individuo->crossover = 0;
	individuo->elite = 0;

	for(int j = 0; j < qtdGenes; ++j)
	{
		individuo->cromossomo[j] = ( rand() % 2 );
	} 
}

//============================================================================================================

void Mutacao(Cromossomo *populacao, int qtdMutacao, int qtdSobreviventes, int tamPopulacao, int qtdGenes)
{
	int indMutacao;
	int indGeneMutacao;

	for(int i = 0; i < qtdMutacao; ++i)
	{
		do
		{
			indMutacao = ( rand() % ( tamPopulacao - 1 ) );
		}
		while( populacao[indMutacao].elite == 1);

		indGeneMutacao = ( rand() % ( qtdGenes - 1 ) );

		populacao[indMutacao].cromossomo[indGeneMutacao] = !populacao[indMutacao].cromossomo[indGeneMutacao];
	}
}

//============================================================================================================

void Crossover(Cromossomo *populacao, int qtdSobreviventes, int tamPopulacao, int qtdGenes, int qtdRenovacao)
{
	int indPai;
	int indMae;
	int indFilho;
	int filhos = ( tamPopulacao - qtdSobreviventes - qtdRenovacao );
	int cont = 0;
	int metadeGenes = ( qtdGenes - ( qtdGenes % 2) ) /2;
	int indSelCrossover = 0;
	
	while(cont < filhos)
	{
		indPai = indSelCrossover++;
		
		if(indSelCrossover >= ( qtdSobreviventes - 1 )) //Volto ao inicio
			indSelCrossover = 0;

		indMae = indSelCrossover++;
		
		//Seleciono um indice para filho desde que nao seja da elite
		do
		{
			indFilho = qtdSobreviventes + ( rand() % ( tamPopulacao - qtdSobreviventes - 1 ) );
		}
		while(populacao[indFilho].elite == 1);
		
		// 50% pai + 50% mae	
		for(int j = 0; j < qtdGenes; ++j)
		{
			if(j < metadeGenes)
				populacao[indFilho].cromossomo[j] = populacao[indPai].cromossomo[j];
			else
				populacao[indFilho].cromossomo[j] = populacao[indMae].cromossomo[j];
		}
		
		++cont;
		
		if(cont < filhos)
		{
			//Seleciono um indice para filho desde que nao seja da elite
			do
			{
				indFilho = qtdSobreviventes + ( rand() % ( tamPopulacao - qtdSobreviventes - 1 ) );
			}
			while(populacao[indFilho].elite == 1);
			
			// 50% mae + 50% pai
			for(int j = 0; j < qtdGenes; ++j)
			{
				if(j < metadeGenes)
					populacao[indFilho].cromossomo[j] = populacao[indMae].cromossomo[j];
				else
					populacao[indFilho].cromossomo[j] = populacao[indPai].cromossomo[j];
			}

			++cont;
		}
	}

	for(int i = 0; i < qtdRenovacao; ++i)
	{
		//Seleciono um indice para filho desde que nao seja da elite
		do
		{
			indFilho = ( rand() % tamPopulacao );
		}
		while(populacao[indFilho].elite == 1);

		CromossomoAleatorio(&populacao[indFilho], qtdGenes);
	}
}

//============================================================================================================

void Selecao(Cromossomo *populacao, int qtdSobreviventes, int tamPopulacao, int qtdElite)
{
	
	float chances;
	float maxChance;

	maxChance = populacao[0].aptidaoRoleta;

	for(int i = 0; i < qtdSobreviventes; ++i)
	{
		if(i < qtdElite)
			populacao[i].elite = 1; //Vai manter com certeza no proxima eliminacao

		chances = ( (float)rand() / RAND_MAX ); //Aleatorio entre 0.0 e 1.0
		chances *= maxChance; 

		//populacao ordenado por aptidao (asc) e consequentemente por aptidaoRoleta (desc)
		for(int j = 0; j < tamPopulacao; ++j)
		{
			if(populacao[j].aptidaoRoleta < chances)
			{
				populacao[j - 1].crossover = 1;//Selecionado para reproduzir no crossover
				break;
			}
		}
	}		
}

//============================================================================================================

long int Fitness(int *cromossomo, long int *genes, int qtdGenes)
{
	long int peso0 = 0;
	long int peso1 = 0;
	long int pesoTotal = 0;
	long int pesoDif = 0;

	for(int i = 0; i < qtdGenes; ++i)
	{
		pesoTotal += genes[i];

		if(cromossomo[i] == 0)
			peso0 += genes[i];
		else
			peso1 += genes[i];
	}

	pesoDif = peso0 - peso1;

	if(pesoDif < 0)
		pesoDif *= -1;

	return ( ( peso0 == 0 || peso1 == 0) ? pesoTotal : pesoDif );
}

//============================================================================================================

void CalculaAptidao(Cromossomo *populacao, int tamPopulacao, long int *genes, int qtdGenes)
{	
	int maxAptidao = 0;
	
	for(int i = 0; i < tamPopulacao; ++i)
	{
		populacao[i].aptidao = Fitness(populacao[i].cromossomo, genes, qtdGenes);
		populacao[i].elite = 0;
		
		if(maxAptidao < populacao[i].aptidao)
			maxAptidao = populacao[i].aptidao;
	}

	float aptidaoAnt = 0;

	for(int i = 0; i < tamPopulacao; ++i)
	{	
		populacao[i].aptidaoRoleta = ( populacao[i].aptidao == 0 ) ? maxAptidao : ( maxAptidao / populacao[i].aptidao ); //Inversamente proporcional a aptidao 
		populacao[i].aptidaoRoleta += aptidaoAnt; //As chances sao o espaco entre as aptidoes
	}
}

//============================================================================================================

int BestFitness(long int *genes, int qtdGenes)
{
	long int peso = 0;
	
	for(int i = 0; i < qtdGenes; ++i)
		peso += genes[i];

	return ( peso % 2 );
}

//============================================================================================================

Cromossomo *PopulacaoInicial(int tamPopulacao, int qtdGenes)
{
	Cromossomo *populacaoInicial = (Cromossomo *)malloc(tamPopulacao * sizeof(Cromossomo));
	
	for(int i = 0; i < tamPopulacao; ++i)
	{
		populacaoInicial[i].cromossomo = (int *)malloc(qtdGenes * sizeof(int));
		
		CromossomoAleatorio(&populacaoInicial[i], qtdGenes);
	}

	return populacaoInicial;
}

//============================================================================================================

int SolucaoEncontrada(Cromossomo *populacao, int tamPopulacao, int bestFitness)
{
	int indSolucao = -1;

	for(int i = 0; i < tamPopulacao; ++i)
	{
		if(populacao[i].aptidao == bestFitness)
		{
			indSolucao = i;
			break;
		}
	}

	return indSolucao;
}

//============================================================================================================

char *CriaArquivo(int qtdLista)
{
	//So cria se nao existir
	int maxLista = 4000000000;
	//char nome[] = "ListaRand";
	char nome[] = "ListaSeq";
	char *filename = (char*)malloc(strlen(nome) * sizeof(char) + 1); 
	strcpy(filename, nome);
	char tamanhoStr[10];
	sprintf(tamanhoStr, "%i", qtdLista);

	filename = (char*)realloc(filename, ( strlen(filename) + strlen(tamanhoStr) )* sizeof(char) + 1); 
	strcat(filename, tamanhoStr);

	char extensao[] = ".txt";
	
	filename = (char*)realloc(filename, ( strlen(filename) + strlen(extensao) )* sizeof(char) + 1); 
	strcat(filename, extensao);

	FILE *file;
	
	if ((file = fopen(filename, "r")))
		fclose(file);	
	else
	{
		time_t segundos;
		struct tm *data_hora_atual;   
		time(&segundos);  
		data_hora_atual = localtime(&segundos);  

		printf("- Gerando arquivo[ %d:%d:%d ]\n",data_hora_atual->tm_hour, data_hora_atual->tm_min, data_hora_atual->tm_sec);
		clock_t inicio = clock();

		ListaSequencial(1, qtdLista, filename);
		//ListaAleatoria(1, maxLista, qtdLista, filename);
		
		time(&segundos);  
		data_hora_atual = localtime(&segundos);  

		clock_t fim = clock();
		
		printf("- Duracao [ %f segundos ]\n", ((float)(fim - inicio))/CLOCKS_PER_SEC);
	}

	return filename;
}

//============================================================================================================

int main (int argc, char **argv)
{
	int qtdGenes = atoi(argv[1]); // 50000 - 100000 - 500000 - 1000000;
	long int *genes;
	int qtdSobreviventes = ( tamPopulacao * percSobreviventes ) / 100;
	int qtdElite = ( tamPopulacao * percElite ) / 100;
	int qtdMutacao = ( tamPopulacao * percMutacao ) / 100;
	int qtdRenovacao = ( tamPopulacao * percRenovacao ) / 100;
	Cromossomo *populacao;
	int indSolucao;
	int *cromSolucao = (int *)malloc(tamPopulacao * sizeof(int));
	long int fitSolucao = 0;
	
	qtdSobreviventes = ( qtdSobreviventes == 0 ) ? 1 : qtdSobreviventes;
	qtdMutacao = ( qtdMutacao == 0 ) ? 1 : qtdMutacao;

	srand(time(NULL)); //Default em C e sempre a mesma seed. Isso modifica a seed a cada iteracao.

	time_t segundos;
	struct tm *data_hora_atual;   
	time(&segundos);  
	data_hora_atual = localtime(&segundos);  

	char *filename = CriaArquivo(qtdGenes); //So cria se nao existir no diretorio
		
	genes = FileToLista(filename, qtdGenes);

	//qsort(genes, qtdGenes, sizeof(genes), ordena);

	printf("- Inicio processamento [ %d:%d:%d ]:\n  - Genes [%i] Populacao [%i] Sobreviventes [%i%%-%i] \n  - Elite [%i%%-%i] Renovacao [%i%%-%i] Mutacao [%i%%-%i] Repeticoes [%i]\n",data_hora_atual->tm_hour, data_hora_atual->tm_min, data_hora_atual->tm_sec, qtdGenes, tamPopulacao, percSobreviventes, qtdSobreviventes, percElite, qtdElite, percRenovacao, qtdRenovacao, percMutacao, qtdMutacao, repeticoes);
	
	clock_t inicio = clock();

	int bestFitness = BestFitness(genes, qtdGenes);

	populacao = PopulacaoInicial(tamPopulacao, qtdGenes);

	//printPopulacao(populacao, tamPopulacao, qtdGenes);
	//printf("-------------------------------------------\n");
		
	int cont = 0;
	int naoMelhorouFit = 0;

	while(cont++ < repeticoes && naoMelhorouFit++ < limSemMelhorarFit)
	{
		qsort(populacao, tamPopulacao, sizeof(Cromossomo), ordenaAptidao);

		CalculaAptidao(populacao, tamPopulacao, genes, qtdGenes);

		if(cont == 1 || populacao[0].aptidao < fitSolucao)
		{
			cromSolucao = populacao[0].cromossomo;
			fitSolucao = populacao[0].aptidao;
			naoMelhorouFit = 0;
		}

		//printPopulacao(populacao, tamPopulacao, qtdGenes);
		//printf("-------------------------------------------\n");
		indSolucao = SolucaoEncontrada(populacao, tamPopulacao, bestFitness);

		if(indSolucao != -1)
			break;
		
		Selecao(populacao, qtdSobreviventes, tamPopulacao, qtdElite);

		qsort(populacao, tamPopulacao, sizeof(Cromossomo), ordenaCrossover);
	
		Crossover(populacao, qtdSobreviventes, tamPopulacao, qtdGenes, qtdRenovacao);

		Mutacao(populacao, qtdMutacao, qtdSobreviventes, tamPopulacao, qtdGenes);
	}

	clock_t fim = clock();
	
	printf("Solucao [ %f segundos ]\n", ((float)(fim - inicio))/CLOCKS_PER_SEC);

	if(indSolucao != -1)
	{
		printf("- Exata apos %i iteracoes:\n", cont - 1);
		PrintSolucao(populacao[indSolucao].cromossomo, genes, qtdGenes);
	}
	else
	{
		printf("- Aproximada apos %i iteracoes [ sem melhorar fitness - %i iteracoes ]:\n", cont - 1, naoMelhorouFit - 1);
		PrintSolucao(cromSolucao, genes, qtdGenes);
	}
	
	free(genes);
	free(populacao->cromossomo);
	free(populacao);	
}
