#include <time.h>

unsigned long int *FileToLista(char *fileName, int tamanhoLista)
{
	FILE *arquivo = fopen(fileName, "r");
	unsigned long int *lista = (unsigned long int*)malloc(tamanhoLista * sizeof(unsigned long int));

	if(arquivo == NULL) 
	{
		fprintf(stderr, "Erro ao abrir o arquivo %s", fileName);
		return NULL;
	}

	int cont = 0;
	
	while( !feof(arquivo) ) 
	{
		if( fscanf(arquivo,"%lu",&lista[cont++]) == 0)
			cont--;
	}
	
	fclose(arquivo);

	return lista;
}

//=================================================================================================================================

void ListaAleatoria(unsigned long int limiteInferior, unsigned long int limiteSuperior, int tamanhoLista, char* nomeArq)
{
	FILE *out = stdout;
	int cont = 0;
	unsigned long int *lista = (unsigned long int*)malloc(tamanhoLista * sizeof(unsigned long int));
	
	srand(time(NULL));

	out = fopen(nomeArq, "w");
	
	while(cont < tamanhoLista)
	{
		unsigned long int aux = limiteInferior + ( rand() % ( limiteSuperior - 1 ) );
		
		if(IndexLista(aux, lista, tamanhoLista) == -1)
			lista[cont++] = aux;
	}

	for(int i = 0; i < tamanhoLista; ++i)
	{
		fprintf(out, "%li ", lista[i] );
	}

	free(lista);

	fclose(out);
}

//=================================================================================================================================

void ListaSequencial(int limiteInferior, int tamanhoLista, char* nomeArq)
{
	FILE *out = stdout;
	
	out = fopen(nomeArq, "w");
	
	for(int i = 1; i <= tamanhoLista; ++i)
		fprintf(out, "%i ", i );
	
	fclose(out);
}
