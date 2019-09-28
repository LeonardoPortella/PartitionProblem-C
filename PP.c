#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include "Util.h"
#include "GeraTeste.h"
#include <mpi.h>
  
int qtdGenes;
unsigned long int *genes;
int repeticoes;
int limSemMelhorFit;
int percMutacao;
int percMutacaoGenes;
int percMutacaoDefault;
int percMutacaoGenesDefault;
int seed;
int geracaoProcMPI = 10;

const int tamPopulacao = 150;
const int qtdNovaGeracao = 50;

//============================================================================================================

typedef struct
{
    int *genotipo;
    unsigned long int aptidao;
    unsigned long int aptidaoRoleta;
	char tipo[10];

} Cromossomo;

//============================================================================================================

typedef struct
{
    Cromossomo *cromossomo;
	Cromossomo melhorCromossomo;
	unsigned long int chancesUltimo;
	int iteracoesSemMelhorFit;
} Populacao;

//============================================================================================================

void ParametrosIniciais(Populacao *populacao)
{
	//Fixo
	repeticoes = 100000;//100000;
	limSemMelhorFit = 500;//500;
	populacao->iteracoesSemMelhorFit = 0;
	percMutacaoGenesDefault = 1;
	percMutacaoDefault = 10;
	percMutacaoGenes = percMutacaoGenesDefault;
	percMutacao = percMutacaoDefault;
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

void PrintSolucao( int *genotipo , unsigned long int *genes , int qtdGenes )
{
    int *solucao0 = ( int * ) malloc( sizeof( int ));
    int tamSol0 = 0;
    unsigned long int peso0 = 0;
    int *solucao1 = ( int * ) malloc( sizeof( int ));;
    int tamSol1 = 0;
    unsigned long int peso1 = 0;
    unsigned long int pesoDif = 0;

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

    if ( peso0 > peso1 )
    	pesoDif = peso0 - peso1;
	else
		pesoDif = peso1 - peso0;

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

	for ( int i = 0 ; i < tamPopulacao ; ++i )
    {
        printf( "- genotipo " );
        PrintLista( populacao->cromossomo[ i ].genotipo , qtdGenes );
        printf( "apt %li " , populacao->cromossomo[ i ].aptidao );
        printf( "aptRoleta %li\n" , populacao->cromossomo[ i ].aptidaoRoleta );
    }
}

//============================================================================================================

unsigned long int Fitness( int *genotipo , unsigned long int *genes , int qtdGenes )
{
    unsigned long int peso0 = 0;
    unsigned long int peso1 = 0;
    unsigned long int pesoDif = 0;

    for ( int i = 0 ; i < qtdGenes ; ++i )
    {
        if ( genotipo[ i ] == 0 )
            peso0 += genes[ i ];
        else
            peso1 += genes[ i ];
    }

	if(peso0 == 0 || peso1 == 0)
		pesoDif = ULONG_MAX;
	else
	{
		if ( peso0 > peso1 )
			pesoDif = peso0 - peso1;
		else
			pesoDif = peso1 - peso0;
	}

	return pesoDif;
}

//============================================================================================================

void cromossomoAleatorio( Cromossomo *individuo , int qtdGenes, unsigned long int *genes, unsigned long int *fitCromossomosNovos, int indAtual, int fitUnico )
{
	int achouIgual;
    individuo->aptidaoRoleta = ULONG_MAX;

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

    unsigned long int maxLista = ULONG_MAX;
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

Populacao *NewPopulacao()
{
	Populacao *populacao = ( Populacao *) malloc( sizeof( Populacao ));

	ParametrosIniciais(populacao);
	populacao->cromossomo = ( Cromossomo * ) malloc( tamPopulacao * sizeof( Cromossomo ));

    populacao->melhorCromossomo.genotipo = ( int * ) malloc( qtdGenes * sizeof( int ));

	for ( int i = 0 ; i < tamPopulacao ; ++i )
        populacao->cromossomo[ i ].genotipo = ( int * ) malloc( qtdGenes * sizeof( int ));

    return populacao;
}

//============================================================================================================

void CromossomoMPI(Cromossomo *cromossomo)
{		
	for(int i = 0; i < qtdGenes; ++i)
		cromossomo->genotipo[ i ] = rand() % 2;

	cromossomo->aptidao = Fitness( cromossomo->genotipo , genes , qtdGenes );

	strcpy(cromossomo->tipo, "Pop ini");
}

//============================================================================================================

void SetCromossomosAleatorios(Populacao *populacao)
{
	printf("Gerando populacao inicial\n");

	int tentativas = 0;

	populacao->chancesUltimo = -1;

	unsigned long int *fitCromossomosNovos = (unsigned long int *)malloc( tamPopulacao * sizeof(unsigned long int));

	for ( int i = 0 ; i < tamPopulacao ; ++i )
    {
		if(tentativas++ <= -1)    	
			cromossomoAleatorio( &populacao->cromossomo[ i ] , qtdGenes, genes, fitCromossomosNovos, i, 1);
		else
			cromossomoAleatorio( &populacao->cromossomo[ i ] , qtdGenes, genes, fitCromossomosNovos, i, 0);
			
		strcpy(populacao->cromossomo[ i ].tipo, "Pop ini");
	}

	free(fitCromossomosNovos);
}

//============================================================================================================

void CalculaAptidao( Populacao *populacao )
{
    qsort( populacao->cromossomo , tamPopulacao, sizeof( Cromossomo ) , ordenaAptidao );

	unsigned long int maxAptidao = 0;
	unsigned long int aptRol = 0;

	for(int i = ( tamPopulacao - 1 ); i >= 0; --i)
	{
        if(populacao->cromossomo[ i ].aptidao != ULONG_MAX)
        {
            maxAptidao = populacao->cromossomo[ i ].aptidao;
            break;
        }
	}

    float aptidaoAnt = 0;

    for ( int i = ( tamPopulacao - 1 ) ; i >= 0 ; --i )
    {
		//Conversão: FaptRol(apt) = max - apt
		//max = 100
		//apt = 0   aptRol = 100
		//apt = 20  aptRol = 80
		//apt = 50  aptRol = 50
		//apt = 79  aptRol = 21
		//apt = 100 aptRol = 0

		if(populacao->cromossomo[ i ].aptidao == ULONG_MAX) //Invalido
			populacao->cromossomo[ i ].aptidaoRoleta = 0;
		else
		{
			aptRol = maxAptidao - populacao->cromossomo[ i ].aptidao;

        	populacao->cromossomo[ i ].aptidaoRoleta = aptRol;
        	populacao->cromossomo[ i ].aptidaoRoleta += aptidaoAnt; //As chances sao o espaco entre as aptidoes

        	aptidaoAnt += aptRol;
		}

		//Armazeno para, durante a Selecao, adicionar as chances do cromossomo com maior aptidao
		if( i == tamPopulacao - 1)
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
		
		memcpy(populacao->melhorCromossomo.genotipo, populacao->cromossomo[ 0 ].genotipo, ( qtdGenes * sizeof(int) ));

		populacao->melhorCromossomo.aptidao = populacao->cromossomo[ 0 ].aptidao;
		//=======================================================================================================================

		strcpy(populacao->melhorCromossomo.tipo, populacao->cromossomo[ 0 ].tipo);

		printf("- Melhorou apos %i iteracoes [ Fitness %li Tipo %s ]\n", populacao->iteracoesSemMelhorFit, populacao->melhorCromossomo.aptidao, populacao->melhorCromossomo.tipo);

		populacao->iteracoesSemMelhorFit = 0;

		percMutacaoGenes = percMutacaoGenesDefault;
		percMutacao = percMutacaoDefault;
	}
	else
	{
		if(iteracoes % 50 == 0)
		{
			++percMutacaoGenes;
			++percMutacao;
		}
	}
}

//============================================================================================================

int *Selecao( Populacao *populacao )
{
	int *paisNovaGeracao = (int *)malloc(qtdNovaGeracao * sizeof(int));
		
	unsigned long int sorteio;
    unsigned long int maxChance;
	
    //populacao ordenado por aptidao (asc) e consequentemente por aptidaoRoleta (desc)
	
    for(int i = 0; i < tamPopulacao; ++i)
	{
		if(populacao->cromossomo[ i ].aptidaoRoleta != 0 ) //Invalido
		{
			maxChance = populacao->cromossomo[ i ].aptidaoRoleta + populacao->chancesUltimo;
			break;
		}
	}

	for ( int i = 0 ; i < qtdNovaGeracao ; ++i )
    {
		sorteio = rand( ) % maxChance;
        paisNovaGeracao[ i ] = 0;

		for ( int j = 0 ; j < tamPopulacao; ++j )
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

void FilhosCrossover(Populacao *populacao, Cromossomo *filho1, Cromossomo *filho2, Cromossomo pai, Cromossomo mae)
{
	//qtdGenes = 1000
	//qtdQuebraGenes1 = 1 + [ 0, 997 ]
	int qtdQuebraGenes1 = 1 + ( rand() % ( qtdGenes - 3 ) ); //0 | 1 2 3 4 5 6 7 8 | 9; 0 1 2 | 3 | 4 5 6 7 8 9
	//qtdQuebraGenes2 = 998 + [ 0, 1 ]
	int qtdQuebraGenes2 = qtdQuebraGenes1 + 1 + ( rand() % ( qtdGenes - qtdQuebraGenes1 - 2) );

	// Percentual de proporcao randomico ( % pai + % mae + % pai)
    for ( int j = 0 ; j < qtdGenes ; ++j )
    {
	    if ( j < qtdQuebraGenes1 )
            filho1->genotipo[ j ] = pai.genotipo[ j ];
        else if( j < qtdQuebraGenes2 )
            filho1->genotipo[ j ] = mae.genotipo[ j ];
        else
            filho1->genotipo[ j ] = pai.genotipo[ j ];
    }

    filho1->aptidao = Fitness( filho1->genotipo , genes , qtdGenes );

	strcpy(filho1->tipo, "Crossover");

	// Percentual de proporcao randomico ( % mae  + % pai + % mae)
    for ( int j = 0 ; j < qtdGenes ; ++j )
    {
	    if ( j < qtdQuebraGenes1 )
            filho2->genotipo[ j ] = pai.genotipo[ j ];
        else if( j < qtdQuebraGenes2 )
            filho2->genotipo[ j ] = mae.genotipo[ j ];
        else
            filho2->genotipo[ j ] = pai.genotipo[ j ];
    }

    filho2->aptidao = Fitness( filho2->genotipo , genes , qtdGenes );

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
    int cont = qtdNovaGeracao;

	for(int i = 0; i < ( qtdNovaGeracao - 1 ); i += 2)
	{
		indPai = paisNovaGeracao[i];
		indMae = paisNovaGeracao[i + 1];

		indFilho1 = tamPopulacao - cont--; //Filhos serao os ultimos na lista de cromossomos
		
		indFilho2 = tamPopulacao - cont--; //Filhos serao os ultimos na lista de cromossomos
		
		//printf("Aqui indPai %i indMae %i indFilho1 %i indFilho2 %i\n", indPai, indMae, indFilho1, indFilho2);

		FilhosCrossover(populacao, &populacao->cromossomo[ indFilho1 ], &populacao->cromossomo[ indFilho2 ], populacao->cromossomo[ indPai ], populacao->cromossomo[ indMae ]);
	}

	free(paisNovaGeracao);
}

//============================================================================================================

void Mutacao( Populacao *populacao )
{
    int indMutacao;
    int qtdGenesMutacao = ( qtdGenes * percMutacaoGenes ) / 100;
	int inicioMutacao;
	int qtdMutacao = ( tamPopulacao * percMutacao ) / 100;
	
	qtdMutacao = ( qtdMutacao == 0 ) ? 1 : qtdMutacao;

	if(qtdGenesMutacao == 0)
        qtdGenesMutacao = 1;
	
	for ( int i = 0 ; i < qtdMutacao ; ++i )
    {
		//Mutacao so ocorrera na nova geracao
        indMutacao = tamPopulacao - ( rand( ) % qtdNovaGeracao ) - 1;

		inicioMutacao = ( rand() % ( qtdGenes - qtdGenesMutacao - 1) );

		for(int j = inicioMutacao; j < ( inicioMutacao + qtdGenesMutacao ); ++j)
			populacao->cromossomo[ indMutacao ].genotipo[ j ] = !populacao->cromossomo[ indMutacao ].genotipo[ j ];

        populacao->cromossomo[ indMutacao ].aptidao = Fitness( populacao->cromossomo[ indMutacao ].genotipo , genes , qtdGenes );

		strcpy(populacao->cromossomo[ indMutacao ].tipo, "Mutacao");
    }
}

//============================================================================================================

void SincronizaPopulacao(Populacao *popResultante, Populacao *popMPI)
{
	//A populacao resultante tera os melhores fitness gerados por cada processo paralelo, 50% de cada um.

	int metade = ( tamPopulacao - ( tamPopulacao % 2 ) ) / 2;

	qsort( popResultante->cromossomo , tamPopulacao , sizeof( Cromossomo ) , ordenaAptidao );

    qsort( popMPI->cromossomo , tamPopulacao , sizeof( Cromossomo ) , ordenaAptidao );

	for(int i = metade; i < tamPopulacao; ++i)
		popResultante->cromossomo[i] = popMPI->cromossomo[i];

	qsort( popResultante->cromossomo , tamPopulacao , sizeof( Cromossomo ) , ordenaAptidao );
}

//============================================================================================================

void ProcessaGeneticoMPI(Populacao *populacao)
{
    int cont = 0;

    while ( cont++ < repeticoes && populacao->iteracoesSemMelhorFit++ < limSemMelhorFit )
    {
		//printf( "CalculaAptidao [ %i ]\n", cont);
		CalculaAptidao( populacao );

		//printf( "Avaliacao [ %i ]\n", cont);
		Avaliacao(populacao, cont);

		if(populacao->melhorCromossomo.aptidao == 0) //Solucao otima (pode nao ser possivel...)
			break;

		//printf( "Selecao e Crossover [ %i ]\n", cont);
		Crossover( populacao, Selecao( populacao ) );

        //printf( "Mutacao [ %i ]\n", cont);
		Mutacao( populacao );
	}
}

//============================================================================================================

int MPI( int argc , char **argv )
{
	int numtasks, pid, rc, dest, source, count, entrada, tag = 1;
	
	MPI_Status status;

	rc = MPI_Init(&argc, &argv);
    
	if(rc != MPI_SUCCESS)
	{
		fprintf(stderr, "Erro ao iniciar o MPI");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	else
	{
		//=====================================================================================
		//Tipo MPI_Cromossomo =================================================================
		//=====================================================================================
		
		int itemsCrom = 4;
		int tamanhosCrom[4] = {1,1,1,10};
		MPI_Datatype tiposCrom[4] = { MPI_INT, MPI_LONG, MPI_LONG, MPI_CHAR };
		MPI_Datatype MPI_Cromossomo;
		MPI_Aint offsetsCrom[4];

		offsetsCrom[0] = offsetof(Cromossomo, genotipo);
		offsetsCrom[1] = offsetof(Cromossomo, aptidao);
		offsetsCrom[2] = offsetof(Cromossomo, aptidaoRoleta);
		offsetsCrom[3] = offsetof(Cromossomo, tipo);

		MPI_Type_create_struct(itemsCrom, tamanhosCrom, offsetsCrom, tiposCrom, &MPI_Cromossomo);
		MPI_Type_commit(&MPI_Cromossomo);

		//=====================================================================================
		//Tipo MPI_Populacao =================================================================
		//=====================================================================================
		/*
	    int itemsPop = 4;
		int tamanhosPop[4] = {1,1,1,1};
		MPI_Datatype tiposPop[4] = { MPI_Cromossomo, MPI_Cromossomo, MPI_LONG, MPI_INT };
		MPI_Datatype MPI_Populacao;
		MPI_Aint offsetsPop[4];

		offsetsPop[0] = offsetof(Populacao, cromossomo);
		offsetsPop[1] = offsetof(Populacao, melhorCromossomo);
		offsetsPop[2] = offsetof(Populacao, chancesUltimo);
		offsetsPop[3] = offsetof(Populacao, iteracoesSemMelhorFit);

		MPI_Type_create_struct(itemsPop, tamanhosPop, offsetsPop, tiposPop, &MPI_Populacao);
		MPI_Type_commit(&MPI_Populacao);
		*/
		//=====================================================================================

		MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);

		if(pid == 0)
		{  
			time_t segundos;
			struct tm *data_hora_atual;
			time( &segundos );
			data_hora_atual = localtime( &segundos );

			clock_t inicio = clock( );

			//=================
			//Populacao Inicial
			//=================
			
			printf("Gerando populacao inicial\n");

			Populacao *populacao = NewPopulacao();

			ParametrosIniciais(populacao);

			int sobras = tamPopulacao % ( numtasks - 1);
			int qtdGerar = ( tamPopulacao - sobras ) / ( numtasks - 1); 

			int aux = qtdGerar;

			for(int i = 1; i < numtasks; ++i)
				MPI_Send(&aux, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
			
			Cromossomo *buffer = ( Cromossomo * )malloc( qtdGerar * sizeof(Cromossomo) );
			
			for(int i = 0; i < qtdGerar; ++i)
				buffer[ i ].genotipo = (int *)malloc(qtdGenes * sizeof(int));
			
			int indPopulacao = 0;

			for(int i = 1; i < numtasks; ++i)
			{
				MPI_Recv(buffer, qtdGerar, MPI_Cromossomo, i, tag, MPI_COMM_WORLD, &status);

				for(int j = 0; j < qtdGerar; ++j)
				{
					populacao->cromossomo[ indPopulacao ] = buffer[ j ];

					populacao->cromossomo[ indPopulacao ].genotipo = (int *)malloc(qtdGenes * sizeof(int));
					MPI_Recv(populacao->cromossomo[ indPopulacao ].genotipo, qtdGenes, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

					++indPopulacao;
				}
			}
			
			free(buffer);

			for(int i = 0; i < sobras; ++i)
			{
				populacao->cromossomo[ indPopulacao ].genotipo = (int *)malloc(qtdGenes * sizeof(int));	
	
				CromossomoMPI(&populacao->cromossomo[ indPopulacao++ ]);
			}

			//=============
			//Processamento
			//=============

			printf
			(
				"- Processamento [ %d:%d:%d ]:\n  - Genes [%i] Populacao [%i] Nova Geracao [%i] \n Mutacao [%i%%] Repeticoes [%i]\n" ,
				    data_hora_atual->tm_hour ,
					data_hora_atual->tm_min ,
					data_hora_atual->tm_sec ,
					qtdGenes , tamPopulacao ,
				    qtdNovaGeracao ,
					percMutacao ,
					repeticoes
			);

			ProcessaGeneticoMPI( populacao );

    		clock_t fim = clock( );

			printf( "Solucao [ %f segundos ]\n" , (( float ) ( fim - inicio )) / CLOCKS_PER_SEC );

			PrintSolucao( populacao->melhorCromossomo.genotipo , genes , qtdGenes );
		
			free( genes );
			free( populacao->cromossomo );
			free( populacao->melhorCromossomo.genotipo );
			free( populacao );
		}
		else if(pid < numtasks)
		{
			//Se todos os processos usarem a mesma seed entao todos gerarao as mesmas populacoes			
			srand(seed + pid);

			//=================
			//Populacao Inicial
			//=================
			
			int qtdGerar;
			
			MPI_Recv(&qtdGerar, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			
			Cromossomo *cromossomo = (Cromossomo *)malloc(qtdGerar * sizeof(Cromossomo));

			for(int i = 0; i < qtdGerar; ++i)
			{
				cromossomo[ i ].genotipo = (int *)malloc(qtdGenes * sizeof(int));	
				CromossomoMPI(&cromossomo[ i ]);
			}
			
			//========================================================================================
			//Obs: Como genotipo e ponteiro dentro de ponteiro, o segundo malloc nao pode ser acessado 
			//fora desta thread e por isso vou precisar enviar separadamente depois
			//========================================================================================
			
			MPI_Send(cromossomo, qtdGerar, MPI_Cromossomo, 0, tag, MPI_COMM_WORLD);

			for(int i = 0; i < qtdGerar; ++i)
				MPI_Send(cromossomo[ i ].genotipo, qtdGenes, MPI_INT, 0, tag, MPI_COMM_WORLD);
			
			//========================================================================================
			
			for(int i = 0; i < qtdGerar; ++i)
				free(cromossomo[ i ].genotipo);
			
			free(cromossomo);

			//=============
			//Processamento
			//=============
		}
		
		MPI_Type_free(&MPI_Cromossomo);
	}
    
	MPI_Finalize();
    }

//============================================================================================================

int main( int argc , char **argv )
{
	if(argc < 3)
	{
		fprintf(stderr, "Informe o tamanho da amostra e a semente\n");
		return -1;
	}
		
	seed = atoi(argv[2]);
	srand( seed );

	qtdGenes = atoi(argv[1]); 
	
	char *filename = CriaArquivo( qtdGenes ); //So cria se nao existir no diretorio

	genes = FileToLista( filename , qtdGenes );
	
	MPI(argc, argv);
}

