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
    int crossover;
    int elite;
    long int aptidaoRoleta;
	int novo;

} Cromossomo;

//============================================================================================================

typedef struct
{
    Cromossomo *cromossomo;
	Cromossomo melhorCromossomo;
	int qtdGenes;
	long int *genes;
	int achouSolucaoOtima;
	long int maxChances;
	int tamPopulacao;
	int repeticoes;
	int limSemMelhorarFit;
	int percSobreviventes;
	int percMutacao;
	int percElite;
	int percRenovacao;
	int qtdSobreviventes;
    int qtdElite;
    int qtdMutacao;
    int qtdRenovacao;
	int qtdFilhos;
	int percEliteMin;
	int percEliteMax;
	int percRenovacaoMin;
	int percRenovacaoMax;
	int percMutacaoMin;
	int percMutacaoMax;
	int rodadasAjusteParam;

} Populacao;

//============================================================================================================

void ParametrosIniciais(Populacao *populacao)
{
	populacao->tamPopulacao = 100;
	populacao->repeticoes = 100000;
	populacao->limSemMelhorarFit = 5000;
	populacao->percSobreviventes = 25;
	populacao->percMutacao = 1;
	populacao->percElite = 50;
	populacao->percRenovacao = 1;
	populacao->rodadasAjusteParam = 2;

	//Barreiras:
	populacao->percEliteMin = 1;
	populacao->percEliteMax = 40;
	populacao->percRenovacaoMin = 1;
	populacao->percRenovacaoMax = 35;
	populacao->percMutacaoMin = 1;
	populacao->percMutacaoMax = 25;
}

//============================================================================================================

void AtualizaParametros(Populacao *populacao)
{
	if(populacao->percElite < populacao->percEliteMin)
		populacao->percElite = populacao->percEliteMin;

	if(populacao->percElite > populacao->percEliteMax)
		populacao->percElite = populacao->percEliteMax;

	if(populacao->percRenovacao < populacao->percRenovacaoMin)
		populacao->percRenovacao = populacao->percRenovacaoMin;

	if(populacao->percRenovacao > populacao->percRenovacaoMax)
		populacao->percRenovacao = populacao->percRenovacaoMax;

	populacao->qtdSobreviventes = ( populacao->tamPopulacao * populacao->percSobreviventes ) / 100;
    populacao->qtdElite = ( populacao->tamPopulacao * populacao->percElite ) / 100;
    populacao->qtdMutacao = ( populacao->tamPopulacao * populacao->percMutacao ) / 100;
    populacao->qtdRenovacao = ( populacao->tamPopulacao * populacao->percRenovacao ) / 100;
    	
	populacao->qtdSobreviventes = ( populacao->qtdSobreviventes == 0 ) ? 1 : populacao->qtdSobreviventes;
    populacao->qtdMutacao = ( populacao->qtdMutacao == 0 ) ? 1 : populacao->qtdMutacao;
	populacao->qtdElite = ( populacao->qtdElite == 0 ) ? 1 : populacao->qtdElite;

	populacao->qtdFilhos = ( populacao->tamPopulacao - populacao->qtdSobreviventes - populacao->qtdRenovacao );
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

int ordenaCrossover( const void *a , const void *b )
{
    Cromossomo *cromossomo0 = ( Cromossomo * ) a;
    Cromossomo *cromossomo1 = ( Cromossomo * ) b;

    return ( cromossomo1->aptidao - cromossomo0->aptidao); //Desc
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
	printf( "Elite %i Crossover %i Novo %i " , populacao->melhorCromossomo.elite, populacao->melhorCromossomo.crossover, populacao->melhorCromossomo.novo );
    printf( "apt %li " , populacao->melhorCromossomo.aptidao );
    printf( "aptRoleta %li\n" , populacao->melhorCromossomo.aptidaoRoleta );
         
	for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
    {
        printf( "- genotipo " );
        PrintLista( populacao->cromossomo[ i ].genotipo , qtdGenes );
        printf( "Elite %i Crossover %i Novo %i " , populacao->cromossomo[ i ].elite, populacao->cromossomo[ i ].crossover, populacao->cromossomo[ i ].novo );
        printf( "apt %li " , populacao->cromossomo[ i ].aptidao );
        printf( "aptRoleta %li\n" , populacao->cromossomo[ i ].aptidaoRoleta );
    }
}

//============================================================================================================

void cromossomoAleatorio( Cromossomo *individuo , int qtdGenes )
{
    individuo->aptidao = -1;
    individuo->aptidaoRoleta = -1;
    individuo->crossover = 0;
    individuo->elite = 0;
	individuo->novo = 1;

    for ( int j = 0 ; j < qtdGenes ; ++j )
    	individuo->genotipo[ j ] = ( rand( ) % 2 );
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

int SolucaoEncontrada( Populacao *populacao , int idealFitness )
{
    int indSolucao = -1;

    for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
    {
        if ( populacao->cromossomo[ i ].aptidao == idealFitness )
        {
            indSolucao = i;
            break;
        }
    }

    return indSolucao;
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
	populacao->cromossomo = ( Cromossomo * ) malloc( populacao->tamPopulacao * sizeof( Cromossomo ));
	
	for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
        populacao->cromossomo[ i ].genotipo = ( int * ) malloc( populacao->qtdGenes * sizeof( int ));

	return populacao;
}

//============================================================================================================

void SetCromossomosAleatorios(Populacao *populacao)
{
	AtualizaParametros(populacao);
	
	populacao->achouSolucaoOtima = -1;
	populacao->maxChances = -1;

	for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
    	cromossomoAleatorio( &populacao->cromossomo[ i ] , populacao->qtdGenes );
}

//============================================================================================================

void CalculaAptidao( Populacao *populacao )
{
    int maxAptidao = 0;
	int aptRol = 0;

	//populacao esta ordenada por aptidao (asc) 

    for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
    {
		if(populacao->cromossomo[ i ].novo == 1) //Calculo o fitness somente se for novo cromossomo
		{
        	populacao->cromossomo[ i ].aptidao = Fitness( populacao->cromossomo[ i ].genotipo , populacao->genes , populacao->qtdGenes );
			populacao->cromossomo[ i ].novo = 0;
        }

		populacao->cromossomo[ i ].elite = 0;

        if ( maxAptidao < populacao->cromossomo[ i ].aptidao )
            maxAptidao = populacao->cromossomo[ i ].aptidao;
    }

    float aptidaoAnt = 0;

    for ( int i = 0 ; i < populacao->tamPopulacao ; ++i )
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
		}

		//Armazeno para, durante a Selecao, adicionar as chances do cromossomo com maior aptidao
		if( i == populacao->tamPopulacao - 1)
			populacao->maxChances = aptRol;
    }
}

//============================================================================================================

void Selecao( Populacao *populacao )
{

    int chances;
    int maxChance;

    //populacao ordenado por aptidao (asc) e consequentemente por aptidaoRoleta (desc)

    maxChance = populacao->cromossomo[ 0 ].aptidaoRoleta;
	
	//Ajuste no criterio da roleta para comtemplar as chances do cromossomo de maior aptidao.	
	//maxChance += populacao->maxChances;

    for ( int i = 0 ; i < populacao->qtdSobreviventes ; ++i )
    {
		if ( i < populacao->qtdElite )
            populacao->cromossomo[ i ].elite = 1; //Vai manter com certeza no proxima eliminacao

        chances = rand( ) % 100; //Aleatorio entre 0 e 100
        chances *= maxChance;

		for ( int j = 0 ; j < populacao->tamPopulacao ; ++j )
        {
			//chance rand [a] = 110; chance rand [b] = 55; chance rand [c] = 20;  chance rand [d] = 50
			//aptRoleta[0] = 100 - [a] -     -     - 
			//aptRoleta[1] = 50  -     - [b] -     - 
			//aptRoleta[2] = 30  -     -     - [c] - [d]
			
            if ( chances > populacao->cromossomo[ j ].aptidaoRoleta || ( j == populacao->tamPopulacao - 1 ))
            {
                populacao->cromossomo[ j ].crossover = 1;//Selecionado para reproduzir no crossover
		        break;
            }
        }
    }
}

//============================================================================================================

void Crossover( Populacao *populacao )
{
    int indPai;
    int indMae;
    int indFilho;
    int cont = 0;
    int metadeGenes = ( populacao->qtdGenes - ( populacao->qtdGenes % 2 )) / 2;
    int indSelCrossover = 0;

    while ( cont < populacao->qtdFilhos )
    {
		indPai = indSelCrossover++;

        if ( indSelCrossover >= ( populacao->qtdSobreviventes - 1 )) //Volto ao inicio
            indSelCrossover = 0;

        indMae = indSelCrossover++;

        //Seleciono um indice para filho desde que nao seja da elite
        do
        {
            indFilho = populacao->qtdSobreviventes + ( rand( ) % ( populacao->tamPopulacao - populacao->qtdSobreviventes - 1 ));
        } while ( populacao->cromossomo[ indFilho ].elite == 1 );

        // 50% pai + 50% mae
        for ( int j = 0 ; j < populacao->qtdGenes ; ++j )
        {
            if ( j < metadeGenes )
                populacao->cromossomo[ indFilho ].genotipo[ j ] = populacao->cromossomo[ indPai ].genotipo[ j ];
            else
                populacao->cromossomo[ indFilho ].genotipo[ j ] = populacao->cromossomo[ indMae ].genotipo[ j ];
        }

        ++cont;

        if ( cont < populacao->qtdFilhos )
        {
            //Seleciono um indice para filho desde que nao seja da elite
            do
            {
                indFilho = populacao->qtdSobreviventes + ( rand( ) % ( populacao->tamPopulacao - populacao->qtdSobreviventes - 1 ));
            } while ( populacao->cromossomo[ indFilho ].elite == 1 );

            // 50% mae + 50% pai
            for ( int j = 0 ; j < populacao->qtdGenes ; ++j )
            {
                if ( j < metadeGenes )
                    populacao->cromossomo[ indFilho ].genotipo[ j ] = populacao->cromossomo[ indMae ].genotipo[ j ];
                else
                    populacao->cromossomo[ indFilho ].genotipo[ j ] = populacao->cromossomo[ indPai ].genotipo[ j ];
            }

            ++cont;
        }
    }

    for ( int i = 0 ; i < populacao->qtdRenovacao ; ++i )
    {
        //Seleciono um indice para filho desde que nao seja da elite
        do
        {
            indFilho = ( rand( ) % populacao->tamPopulacao );
        } while ( populacao->cromossomo[ indFilho ].elite == 1 );

        cromossomoAleatorio( &populacao->cromossomo[ indFilho ] , populacao->qtdGenes );
    }
}

//============================================================================================================

void Mutacao( Populacao *populacao )
{
    int indMutacao;
    int indGeneMutacao;

    for ( int i = 0 ; i < populacao->qtdMutacao ; ++i )
    {
        do
        {
            indMutacao = ( rand( ) % ( populacao->tamPopulacao - 1 ));
        } while ( populacao->cromossomo[ indMutacao ].elite == 1 );

        indGeneMutacao = ( rand( ) % ( populacao->qtdGenes - 1 ));

        populacao->cromossomo[ indMutacao ].genotipo[ indGeneMutacao ] = !populacao->cromossomo[ indMutacao ].genotipo[ indGeneMutacao ];
    }
}

//============================================================================================================

void SincronizaPopulacaoMPI(Populacao *popResultante, Populacao *popMPI)
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

void ProcessaGeneticoMPI(Populacao *populacaoMPI)
{
    int cont = 0;
    int naoMelhorouFit = 0;
	int indSolucao = -1;

    while ( cont++ < populacaoMPI->repeticoes && naoMelhorouFit++ < populacaoMPI->limSemMelhorarFit )
    {
		//printPopulacao(populacaoMPI, qtdGenes);
    	//printf("-------------------------------------------\n");
		//sleep(4);
    	//printf( "CalculaAptidao\n");
		CalculaAptidao( populacaoMPI );

        qsort( populacaoMPI->cromossomo , populacaoMPI->tamPopulacao , sizeof( Cromossomo ) , ordenaAptidao );

        if ( cont == 1 || populacaoMPI->cromossomo[ 0 ].aptidao < populacaoMPI->melhorCromossomo.aptidao )
        {
            populacaoMPI->melhorCromossomo = populacaoMPI->cromossomo[ 0 ];
            naoMelhorouFit = 0;
        
			//printf("MELHOROU\n");
		
			if(populacaoMPI->melhorCromossomo.aptidao == 0) //Solucao otima
				break;

			if(cont % populacaoMPI->rodadasAjusteParam == 0)
			{
				populacaoMPI->percElite++;
				populacaoMPI->percMutacao--;
				populacaoMPI->percRenovacao--;
			}
		}
		else
		{
			if(cont % populacaoMPI->rodadasAjusteParam == 0)
			{
				populacaoMPI->percElite--;
				populacaoMPI->percMutacao++;
				populacaoMPI->percRenovacao++;
			}
		}
		
		AtualizaParametros(populacaoMPI);
			       
		//printf( "Selecao\n");
		Selecao( populacaoMPI );

        qsort( populacaoMPI->cromossomo , populacaoMPI->tamPopulacao , sizeof( Cromossomo ) , ordenaCrossover );

		//printf( "Crossover\n");
		Crossover( populacaoMPI );

		//printf( "Mutacao\n");
		Mutacao( populacaoMPI );
   }

}

//============================================================================================================

void ProcessaGenetico(int qtdGenes)
{
	srand( time(NULL)); //Default em C e sempre a mesma seed. Isso modifica a seed a cada iteracao.

    time_t segundos;
    struct tm *data_hora_atual;
    time( &segundos );
    data_hora_atual = localtime( &segundos );

    clock_t inicio = clock( );

	Populacao *populacao = NewPopulacao(NULL, qtdGenes);
	SetCromossomosAleatorios( populacao );
	
	printf
	( 
		"- Processamento [ %d:%d:%d ]:\n  - Genes [%i] Populacao [%i] Sobreviventes [%i%%-%i] \n  - Elite [%i%%-%i] Renovacao [%i%%-%i] Mutacao [%i%%-%i] Repeticoes [%i]\n" ,
            data_hora_atual->tm_hour , 
			data_hora_atual->tm_min , 
			data_hora_atual->tm_sec , 
			qtdGenes , populacao->tamPopulacao ,
            populacao->percSobreviventes , 
			populacao->qtdSobreviventes , 
			populacao->percElite , 
			populacao->qtdElite , 
			populacao->percRenovacao , 
			populacao->qtdRenovacao , 
			populacao->percMutacao ,
            populacao->qtdMutacao , 
			populacao->repeticoes 
	);

	ProcessaGeneticoMPI(populacao);

    clock_t fim = clock( );

    printf( "Solucao [ %f segundos ]\n" , (( float ) ( fim - inicio )) / CLOCKS_PER_SEC );

    PrintSolucao( populacao->melhorCromossomo.genotipo , populacao->genes , qtdGenes );
    
	free( populacao->genes );	
    free( populacao->cromossomo );
    free( populacao );
}

//============================================================================================================

int main( int argc , char **argv )
{
	if(argc == 1)
	{
		fprintf(stderr, "Informe o tamanho da amostra\n");
		return 1;
	}

	ProcessaGenetico(atoi(argv[1]));
}
