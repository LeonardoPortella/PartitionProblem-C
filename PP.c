#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include "Util.h"
#include "GeraTeste.h"
//#include <mpi.h>

const int instanciasMPI = 10;

//============================================================================================================

typedef struct
{
    int *genotipo;
    long int aptidao;
    long int aptidaoRoleta;
	char tipo[10];

} Cromossomo;

//============================================================================================================

typedef struct
{
    Cromossomo *cromossomo;
	Cromossomo melhorCromossomo;
	int qtdGenes;
	long int *genes;
	long int chancesUltimo;
	int tamPopulacao;
	int repeticoes;
	int limSemMelhorFit;
	int percNovaGeracao;
	int percMutacao;
	int percRenovacao;
	int qtdNovaGeracao;
    int qtdMutacao;
    int qtdRenovacao;
	int qtdFilhos;
	int percRenovacaoMin;
	int percRenovacaoMax;
	int percMutacaoMin;
	int percMutacaoMax;
	int percMutacaoGenes;
	int rodadasAjusteParam;
	int percMutacaoDefault;
	int percRenovacaoDefault;
	int iteracoesSemMelhorFit;

} Populacao;

//============================================================================================================

void ParametrosIniciais(Populacao *populacao)
{
	//Fixo
	populacao->tamPopulacao = 100;
	populacao->repeticoes = 100000;//100000;
	populacao->limSemMelhorFit = 500;//500;
	populacao->rodadasAjusteParam = 1;
	populacao->percNovaGeracao = 50;
	populacao->iteracoesSemMelhorFit = 0;
	populacao->percMutacaoGenes = 1;

	//Default
	populacao->percRenovacaoDefault = 0;
	populacao->percMutacaoDefault = 1;
	populacao->percMutacao = populacao->percMutacaoDefault;
	populacao->percRenovacao = populacao->rodadasAjusteParam;

	//Barreiras:
	populacao->percRenovacaoMin = populacao->percRenovacaoDefault;
	populacao->percRenovacaoMax = 0;
	populacao->percMutacaoMin = populacao->percMutacaoDefault;
	populacao->percMutacaoMax = 10;

	populacao->qtdNovaGeracao = ( populacao->tamPopulacao * populacao->percNovaGeracao ) / 100;
	populacao->qtdNovaGeracao = ( populacao->qtdNovaGeracao == 0 ) ? 1 : populacao->qtdNovaGeracao;
}

//============================================================================================================

void AtualizaParametros(Populacao *populacao)
{
	if(populacao->percRenovacao < populacao->percRenovacaoMin)
		populacao->percRenovacao = populacao->percRenovacaoMin;

	if(populacao->percRenovacao > populacao->percRenovacaoMax)
		populacao->percRenovacao = populacao->percRenovacaoMax;

    populacao->qtdMutacao = ( populacao->tamPopulacao * populacao->percMutacao ) / 100;
    populacao->qtdRenovacao = ( populacao->tamPopulacao * populacao->percRenovacao ) / 100;

    populacao->qtdMutacao = ( populacao->qtdMutacao == 0 ) ? 1 : populacao->qtdMutacao;

	populacao->qtdFilhos = ( populacao->tamPopulacao - populacao->qtdNovaGeracao - populacao->qtdRenovacao );

	if(populacao->limSemMelhorFit > populacao->qtdGenes)
		populacao->limSemMelhorFit = populacao->qtdGenes;
}

//============================================================================================================

int ordena( const void *a , const void *b )
{
    return ( a - b );
}

//============================================================================================================

int ordenaAptidao( const void *a , const void *b )
{
    Cromossomo *cromossomo0 = ( Cromossomo * ) a;
    Cromossomo *cromossomo1 = ( Cromossomo * ) b;

    return ( cromossomo0->aptidao - cromossomo1->aptidao );
}

//============================================================================================================

void PrintSolucao( int *genotipo , long int *genes , int qtdGenes )
{
    int *solucao0 = ( int * ) malloc( sizeof( int ));
    int tamSol0 = 0;
    long int peso0 = 0;
    int *solucao1 = ( int * ) malloc( sizeof( int ));;
    int tamSol1 = 0;
    long int peso1 = 0;
    long int pesoDif = 0;

    for ( int i = 0 ; i < qtdGenes ; ++i )
    {
        if ( genotipo[ i ] == 0 )
        {
            solucao0 = ( int * ) realloc( solucao0 , ++tamSol0 * sizeof( int ));
            solucao0[ tamSol0 - 1 ] = genes[ i ];
            peso0 += genes[ i ];
        } else
        {
            solucao1 = ( int * ) realloc( solucao1 , ++tamSol1 * sizeof( int ));
            solucao1[ tamSol1 - 1 ] = genes[ i ];
            peso1 += genes[ i ];
        }
    }

    printf( "Peso[0] %li " , peso0 );

    if ( qtdGenes <= 15 )
    {
		printf( "Sol[0] ");
	    PrintLista( solucao0 , tamSol0 );
	}

    printf( "Peso[1] %li " , peso1 );

    if ( qtdGenes <= 15 )
	{
    	printf( "Sol[1] ");
	    PrintLista( solucao1 , tamSol1 );
	}

    pesoDif = peso0 - peso1;

    if ( pesoDif < 0 )
        pesoDif *= -1;

    printf( "Diferenca Peso %li\n" , pesoDif );

    free( solucao0 );
    free( solucao1 );
}

//============================================================================================================

void printPopulacao( Populacao *populacao , int qtdGenes )
{
    printf( "- Melhor genotipo " );
    PrintLista( populacao->melhorCromossomo.genotipo , qtdGenes );
    printf( "apt %li " , populacao->melhorCromossomo.aptidao );
    printf( "aptRoleta %li\n" , populacao->melhorCromossomo.aptidaoRoleta );

	for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
    {
        printf( "- genotipo " );
        PrintLista( populacao->cromossomo[ i ].genotipo , qtdGenes );
        printf( "apt %li " , populacao->cromossomo[ i ].aptidao );
        printf( "aptRoleta %li\n" , populacao->cromossomo[ i ].aptidaoRoleta );
    }
}

//============================================================================================================

long int Fitness( int *genotipo , long int *genes , int qtdGenes )
{
    long int peso0 = 0;
    long int peso1 = 0;
    long int pesoDif = 0;

    for ( int i = 0 ; i < qtdGenes ; ++i )
    {
        if ( genotipo[ i ] == 0 )
            peso0 += genes[ i ];
        else
            peso1 += genes[ i ];
    }

	if(peso0 == 0 || peso1 == 0)
		pesoDif = INT_MAX;
	else
	{
    	pesoDif = peso0 - peso1;

		if ( pesoDif < 0 )
		    pesoDif *= -1;
	}

	return pesoDif;
}

//============================================================================================================

void cromossomoAleatorio( Cromossomo *individuo , int qtdGenes, long int *genes, long int *fitCromossomosNovos, int indAtual, int fitUnico )
{
	int achouIgual;
    individuo->aptidaoRoleta = INT_MAX;

	do
	{
		for ( int j = 0 ; j < qtdGenes ; ++j )
		{
			individuo->genotipo[ j ] = ( rand( ) % 2 );

            individuo->aptidao = Fitness( individuo->genotipo , genes , qtdGenes );
		}

		if(fitUnico == 0)
			break;

		achouIgual = 0; //Remove duplicados e simetricos na geracao aleatoria

		for( int i = 0; i < indAtual; ++i)
		{
			if( fitCromossomosNovos[i] == individuo->aptidao )
			{
				achouIgual = 1;
				break;
			}
		}

	} while ( achouIgual == 1 );
}

//============================================================================================================

char *CriaArquivo( int qtdLista )
{
    //So cria se nao existir

    int maxLista = INT_MAX;
    //char nome[] = "ListaRand";
    char nome[] = "ListaSeq";
    char *filename = ( char * ) malloc( strlen( nome ) * sizeof( char ) + 1 );
    strcpy( filename , nome );
    char tamanhoStr[10];
    sprintf( tamanhoStr , "%i" , qtdLista );

    filename = ( char * ) realloc( filename , ( strlen( filename ) + strlen( tamanhoStr )) * sizeof( char ) + 1 );
    strcat( filename , tamanhoStr );

    char extensao[] = ".txt";

    filename = ( char * ) realloc( filename , ( strlen( filename ) + strlen( extensao )) * sizeof( char ) + 1 );
    strcat( filename , extensao );

    FILE *file;

    if (( file = fopen( filename , "r" )))
        fclose( file );
    else
    {
        time_t segundos;
        struct tm *data_hora_atual;
        time( &segundos );
        data_hora_atual = localtime( &segundos );

        printf( "- Gerando arquivo[ %d:%d:%d ]\n" , data_hora_atual->tm_hour , data_hora_atual->tm_min ,
                data_hora_atual->tm_sec );
        clock_t inicio = clock( );

        ListaSequencial( 1 , qtdLista , filename );
        //ListaAleatoria(1, maxLista, qtdLista, filename);

        time( &segundos );
        data_hora_atual = localtime( &segundos );

        clock_t fim = clock( );

        printf( "- Duracao [ %f segundos ]\n" , (( float ) ( fim - inicio )) / CLOCKS_PER_SEC );
    }

    return filename;
}

//============================================================================================================

int Alocou2Vezes(Populacao *populacao)
{

    int **memoria = (int **)malloc( ( populacao->tamPopulacao ) * sizeof(int *) );

    for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
        memoria[ i ] = populacao->cromossomo[ i ].genotipo;

	for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
    {
        //printf("memoria [ %i ] %p\n",i, memoria[ i ]);

        if(memoria[ i ] == populacao->melhorCromossomo.genotipo)
		{
            printf("Alocou 2 vezes [ %p ] indice [ %i ]!!!\n", memoria[ i ], i);
			free(memoria);
			return 1;
		}
    }

	free(memoria);

	return 0;
}

//============================================================================================================

Populacao *NewPopulacao(long int *genes, int qtdGenes)
{
	Populacao *populacao = ( Populacao *) malloc( sizeof( Populacao ));

	if(genes == NULL)
	{
		char *filename = CriaArquivo( qtdGenes ); //So cria se nao existir no diretorio

    	populacao->genes = FileToLista( filename , qtdGenes );
    }
	else
    	populacao->genes = genes;

	populacao->qtdGenes = qtdGenes;

	ParametrosIniciais(populacao);
	populacao->cromossomo = ( Cromossomo * ) malloc( ( populacao->tamPopulacao + populacao->qtdNovaGeracao ) * sizeof( Cromossomo ));

    populacao->melhorCromossomo.genotipo = ( int * ) malloc( populacao->qtdGenes * sizeof( int ));

	for ( int i = 0 ; i < ( populacao->tamPopulacao + populacao->qtdNovaGeracao ) ; ++i )
        populacao->cromossomo[ i ].genotipo = ( int * ) malloc( populacao->qtdGenes * sizeof( int ));

    return populacao;
}

//============================================================================================================

Populacao *NewPopulacaoMPI(long int *genes, int qtdGenes)
{
	Populacao *populacao = ( Populacao *) malloc( sizeof( Populacao ));

	if(genes == NULL)
	{
		char *filename = CriaArquivo( qtdGenes ); //So cria se nao existir no diretorio

    	populacao->genes = FileToLista( filename , qtdGenes );
    }
	else
    	populacao->genes = genes;

	populacao->qtdGenes = qtdGenes;

	ParametrosIniciais(populacao);
	populacao->cromossomo = ( Cromossomo * ) malloc( ( populacao->tamPopulacao + populacao->qtdNovaGeracao ) * sizeof( Cromossomo ));

    populacao->melhorCromossomo.genotipo = ( int * ) malloc( populacao->qtdGenes * sizeof( int ));

	for ( int i = 0 ; i < ( populacao->tamPopulacao + populacao->qtdNovaGeracao ) ; ++i )
        populacao->cromossomo[ i ].genotipo = ( int * ) malloc( populacao->qtdGenes * sizeof( int ));

    return populacao;
}

//============================================================================================================

void SetCromossomosAleatorios(Populacao *populacao)
{
	int tentativas = 0;

	AtualizaParametros(populacao);

	populacao->chancesUltimo = -1;

	long int *fitCromossomosNovos = (long int *)malloc(( populacao->tamPopulacao + populacao->qtdNovaGeracao ) * sizeof(long int));

	for ( int i = 0 ; i < ( populacao->tamPopulacao + populacao->qtdNovaGeracao ) ; ++i )
    {
		if(tentativas++ <= 1)    	
			cromossomoAleatorio( &populacao->cromossomo[ i ] , populacao->qtdGenes, populacao->genes, fitCromossomosNovos, i, 1);
		else
			cromossomoAleatorio( &populacao->cromossomo[ i ] , populacao->qtdGenes, populacao->genes, fitCromossomosNovos, i, 0);
			
		strcpy(populacao->cromossomo[ i ].tipo, "Pop ini");
	}

	free(fitCromossomosNovos);
}

//============================================================================================================

void CalculaAptidao( Populacao *populacao )
{
    qsort( populacao->cromossomo , ( populacao->tamPopulacao + populacao->qtdNovaGeracao ) , sizeof( Cromossomo ) , ordenaAptidao );

	long int maxAptidao = 0;
	long int aptRol = 0;

	for(int i = ( populacao->tamPopulacao + populacao->qtdNovaGeracao - 1 ); i >= 0; --i)
	{
        if(populacao->cromossomo[ i ].aptidao != INT_MAX)
        {
            maxAptidao = populacao->cromossomo[ i ].aptidao;
            break;
        }
	}

    float aptidaoAnt = 0;

    for ( int i = ( populacao->tamPopulacao + populacao->qtdNovaGeracao - 1 ) ; i >= 0 ; --i )
    {
		//Conversão: FaptRol(apt) = max - apt
		//max = 100
		//apt = 0   aptRol = 100
		//apt = 20  aptRol = 80
		//apt = 50  aptRol = 50
		//apt = 79  aptRol = 21
		//apt = 100 aptRol = 0

		if(populacao->cromossomo[ i ].aptidao == INT_MAX) //Invalido
			populacao->cromossomo[ i ].aptidaoRoleta = 0;
		else
		{
			aptRol = maxAptidao - populacao->cromossomo[ i ].aptidao;

        	populacao->cromossomo[ i ].aptidaoRoleta = aptRol;
        	populacao->cromossomo[ i ].aptidaoRoleta += aptidaoAnt; //As chances sao o espaco entre as aptidoes

        	aptidaoAnt += aptRol;
		}

		//Armazeno para, durante a Selecao, adicionar as chances do cromossomo com maior aptidao
		if( i == populacao->tamPopulacao - 1)
			populacao->chancesUltimo = aptidaoAnt;
    }
}

//============================================================================================================

void Avaliacao( Populacao *populacao, int iteracoes)
{
	if ( iteracoes == 1 || populacao->cromossomo[ 0 ].aptidao < populacao->melhorCromossomo.aptidao )
    {
    	//=======================================================================================================================
		//=== Motivo de varias horas de debug: Copiar o struct Cromossomo, que tem um ponteiro de int em sua estrutura, =========
		//=== faz com que na copia, o ponteiro de int tenha o mesmo endereco do ponteiro da copia. ==============================
		//=== Precisa copiar o conteudo do ponteiro interno apenas ==============================================================
		//=== Para copiar estrutura utilizar memcpy =============================================================================
		//=======================================================================================================================
		memcpy(populacao->melhorCromossomo.genotipo, populacao->cromossomo[ 0 ].genotipo, ( populacao->qtdGenes * sizeof(int) ));

		populacao->melhorCromossomo.aptidao = populacao->cromossomo[ 0 ].aptidao;
		//=======================================================================================================================

		strcpy(populacao->melhorCromossomo.tipo, populacao->cromossomo[ 0 ].tipo);

		printf("- Melhorou apos %i iteracoes [ Fitness %li Tipo %s ]\n", populacao->iteracoesSemMelhorFit, populacao->melhorCromossomo.aptidao, populacao->melhorCromossomo.tipo);

		populacao->iteracoesSemMelhorFit = 0;

		populacao->percMutacao = populacao->percMutacaoDefault;
		populacao->percRenovacao = populacao->percRenovacaoDefault;
	}
	else if(iteracoes % populacao->rodadasAjusteParam == 0)
	{
		populacao->percMutacao++;
		populacao->percRenovacao++;

		AtualizaParametros(populacao);
	}
}

//============================================================================================================

int *Selecao( Populacao *populacao )
{
	long int sorteio;
    long int maxChance;
	int *paisNovaGeracao = (int *)malloc(populacao->qtdNovaGeracao * sizeof(int));

    //populacao ordenado por aptidao (asc) e consequentemente por aptidaoRoleta (desc)

    for(int i = 0; i < populacao->tamPopulacao; ++i)
	{
		if(populacao->cromossomo[ i ].aptidaoRoleta != 0 ) //Invalido
		{
			maxChance = populacao->cromossomo[ i ].aptidaoRoleta + populacao->chancesUltimo;
			break;
		}
	}

	for ( int i = 0 ; i < populacao->qtdNovaGeracao ; ++i )
    {
		sorteio = rand( ) % maxChance;
        paisNovaGeracao[ i ] = 0;

		for ( int j = 0 ; j < ( populacao->tamPopulacao + populacao->qtdNovaGeracao ); ++j )
        {
			//sorteio rand [a] = 110; sorteio rand [b] = 55; sorteio rand [c] = 20;  sorteio rand [d] = 50
			//aptRoleta[0] = 100 - [a] -     -     -
			//aptRoleta[1] = 50  -     - [b] -     -
			//aptRoleta[2] = 30  -     -     - [c] - [d]

            if ( sorteio < populacao->cromossomo[ j ].aptidaoRoleta )
                paisNovaGeracao[ i ] = j;//Selecionado para reproduzir no crossover
            else
                break;
        }
    }

	return paisNovaGeracao;
}

//============================================================================================================

void Renovacao( Populacao *populacao)
{
    int indFilho;
    long int *fitCromossomosNovos = (long int *)malloc(populacao->qtdRenovacao * sizeof(long int));

	for ( int i = 0 ; i < populacao->qtdRenovacao ; ++i )
    {
        indFilho = populacao->tamPopulacao + ( rand( ) % ( populacao->qtdNovaGeracao - 1 ));

        cromossomoAleatorio( &populacao->cromossomo[ indFilho ] , populacao->qtdGenes, populacao->genes, fitCromossomosNovos, i, 1 );

		strcpy(populacao->cromossomo[ indFilho ].tipo, "Renovacao");
    }

	free(fitCromossomosNovos);
}

//============================================================================================================

void FilhosCrossover(Populacao *populacao, Cromossomo *filho1, Cromossomo *filho2, Cromossomo pai, Cromossomo mae)
{
	//populacao->qtdGenes = 1000
	//qtdQuebraGenes1 = 1 + [ 0, 997 ]
	int qtdQuebraGenes1 = 1 + ( rand() % ( populacao->qtdGenes - 3 ) ); //0 | 1 2 3 4 5 6 7 8 | 9; 0 1 2 | 3 | 4 5 6 7 8 9
	//qtdQuebraGenes2 = 998 + [ 0, 1 ]
	int qtdQuebraGenes2 = qtdQuebraGenes1 + 1 + ( rand() % ( populacao->qtdGenes - qtdQuebraGenes1 - 2) );

	// Percentual de proporcao randomico ( % pai + % mae + % pai)
    for ( int j = 0 ; j < populacao->qtdGenes ; ++j )
    {
	    if ( j < qtdQuebraGenes1 )
            filho1->genotipo[ j ] = pai.genotipo[ j ];
        else if( j < qtdQuebraGenes2 )
            filho1->genotipo[ j ] = mae.genotipo[ j ];
        else
            filho1->genotipo[ j ] = pai.genotipo[ j ];
    }

    filho1->aptidao = Fitness( filho1->genotipo , populacao->genes , populacao->qtdGenes );

	strcpy(filho1->tipo, "Crossover");

	// Percentual de proporcao randomico ( % mae  + % pai + % mae)
    for ( int j = 0 ; j < populacao->qtdGenes ; ++j )
    {
	    if ( j < qtdQuebraGenes1 )
            filho2->genotipo[ j ] = pai.genotipo[ j ];
        else if( j < qtdQuebraGenes2 )
            filho2->genotipo[ j ] = mae.genotipo[ j ];
        else
            filho2->genotipo[ j ] = pai.genotipo[ j ];
    }

    filho2->aptidao = Fitness( filho2->genotipo , populacao->genes , populacao->qtdGenes );

	strcpy(filho2->tipo, "Crossover");
}

//============================================================================================================

void Crossover( Populacao *populacao, int *paisNovaGeracao )
{
	//Crossover de 2 pontos aleatorios

	int indPai;
    int indMae;
    int indFilho1;
    int indFilho2;
    int cont = 0;

	for(int i = 0; i < ( populacao->qtdNovaGeracao - 1 ); i += 2)
	{
		indPai = paisNovaGeracao[i];
		indMae = paisNovaGeracao[i + 1];

		indFilho1 = cont++ + populacao->tamPopulacao; //Filhos serao os ultimos na lista de cromossomos
		populacao->cromossomo[ indFilho1 ];

		indFilho2 = cont++ + populacao->tamPopulacao; //Filhos serao os ultimos na lista de cromossomos
		populacao->cromossomo[ indFilho2 ];

		FilhosCrossover(populacao, &populacao->cromossomo[ indFilho1 ], &populacao->cromossomo[ indFilho2 ], populacao->cromossomo[ indPai ], populacao->cromossomo[ indMae ]);
	}

	free(paisNovaGeracao);
}

//============================================================================================================

void Mutacao( Populacao *populacao )
{
    int indMutacao;
    int qtdGenesMutacao = ( populacao->qtdGenes * populacao->percMutacaoGenes ) / 100;
	int inicioMutacao;

	if(qtdGenesMutacao == 0)
        qtdGenesMutacao = 1;

    for ( int i = 0 ; i < populacao->qtdMutacao ; ++i )
    {
		//Mutacao so ocorrera na nova geracao
        indMutacao = populacao->tamPopulacao + ( rand( ) % ( populacao->qtdNovaGeracao - 1 ));

		inicioMutacao = ( rand() % ( populacao->qtdGenes - qtdGenesMutacao - 1 ) );

    	for(int j = inicioMutacao; j < ( inicioMutacao + qtdGenesMutacao ); ++j)
			populacao->cromossomo[ indMutacao ].genotipo[ j ] = !populacao->cromossomo[ indMutacao ].genotipo[ j ];

        populacao->cromossomo[ indMutacao ].aptidao = Fitness( populacao->cromossomo[ indMutacao ].genotipo , populacao->genes , populacao->qtdGenes );

		strcpy(populacao->cromossomo[ indMutacao ].tipo, "Mutacao");
    }
}

//============================================================================================================

void SincronizaPopulacao(Populacao *popResultante, Populacao *popMPI)
{
	//A populacao resultante tera os melhores fitness gerados por cada processo paralelo, 50% de cada um.

	int metade = ( popResultante->tamPopulacao - popResultante->tamPopulacao % 2 ) / 2;

	qsort( popResultante->cromossomo , popResultante->tamPopulacao , sizeof( Cromossomo ) , ordenaAptidao );

    qsort( popMPI->cromossomo , popMPI->tamPopulacao , sizeof( Cromossomo ) , ordenaAptidao );

	for(int i = metade; i < popResultante->tamPopulacao; ++i)
		popResultante->cromossomo[i] = popMPI->cromossomo[i];

	qsort( popResultante->cromossomo , popResultante->tamPopulacao , sizeof( Cromossomo ) , ordenaAptidao );
}

//============================================================================================================

void ProcessaGeneticoMPI(Populacao *populacao)
{
    int cont = 0;

    while ( cont++ < populacao->repeticoes && populacao->iteracoesSemMelhorFit++ < populacao->limSemMelhorFit )
    {
		//printf( "CalculaAptidao [ %i ]\n", cont);
		CalculaAptidao( populacao );

        //printf( "Avaliacao [ %i ]\n", cont);
		Avaliacao(populacao, cont);

		if(populacao->melhorCromossomo.aptidao == 0) //Solucao otima (pode nao ser possivel...)
			break;

		//printf( "Selecao/Crossover [ %i ]\n", cont);
		Crossover( populacao, Selecao( populacao ) );

		//printf( "Renovacao [ %i ]\n", cont);
		//Renovacao( populacao );

        //printf( "Mutacao [ %i ]\n", cont);
		Mutacao( populacao );
	}
}

//============================================================================================================

void ProcessaGenetico(int qtdGenes, int seed)
{
	time_t segundos;
    struct tm *data_hora_atual;
    time( &segundos );
    data_hora_atual = localtime( &segundos );

    clock_t inicio = clock( );

	printf("Gerando populacao inicial\n");

	Populacao *populacao = NewPopulacao(NULL, qtdGenes);
	SetCromossomosAleatorios( populacao );

	printf
	(
		"- Processamento [ %d:%d:%d ]:\n  - Genes [%i] Populacao [%i] Nova Geracao [%i%%-%i] \n Renovacao [%i%%-%i] Mutacao [%i%%-%i] Repeticoes [%i]\n" ,
            data_hora_atual->tm_hour ,
			data_hora_atual->tm_min ,
			data_hora_atual->tm_sec ,
			qtdGenes , populacao->tamPopulacao ,
            populacao->percNovaGeracao ,
			populacao->qtdNovaGeracao ,
			populacao->percRenovacao ,
			populacao->qtdRenovacao ,
			populacao->percMutacao ,
            populacao->qtdMutacao ,
			populacao->repeticoes
	);

	ProcessaGeneticoMPI( populacao );

    clock_t fim = clock( );

    printf( "Solucao [ %f segundos ]\n" , (( float ) ( fim - inicio )) / CLOCKS_PER_SEC );

    PrintSolucao( populacao->melhorCromossomo.genotipo , populacao->genes , qtdGenes );

	free( populacao->genes );

    for(int i = 0; i < populacao->tamPopulacao; ++i)
        free( populacao->cromossomo[ i ].genotipo );

    free( populacao->cromossomo );
    free( populacao->melhorCromossomo.genotipo );
    free( populacao );
}

//============================================================================================================

int main( int argc , char **argv )
{
	if(argc < 3)
	{
		fprintf(stderr, "Informe o tamanho da amostra e a semente\n");
		return 1;
	}

	srand( atoi(argv[2]) );

	ProcessaGenetico(atoi(argv[1]), 5);
}

