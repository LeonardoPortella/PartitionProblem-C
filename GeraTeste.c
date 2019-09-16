long int *FileToLista(char *fileName, int tamanhoLista)
{
	FILE *arquivo = fopen(fileName, "r");
	long int *lista = (long int*)malloc(tamanhoLista * sizeof(long int));

	if(arquivo == NULL) 
	{
		fprintf(stderr, "Erro ao abrir o arquivo %s", fileName);
		return NULL;
	}

	int cont = 0;
	
	while( !feof(arquivo) ) 
	{
		fscanf(arquivo,"%li",&lista[cont++]);
	}

	fclose(arquivo);

	return lista;
}

void ListaAleatoria(int limiteInferior, int limiteSuperior, int tamanhoLista, char* nomeArq)
{
	FILE *out = stdout;
	int cont = 0;
	int *lista = (int*)malloc(tamanhoLista * sizeof(int));
	
	out = fopen(nomeArq, "w");
	
	while(cont < tamanhoLista)
	{
		int aux = limiteInferior + ( rand() % ( limiteSuperior - 1 ) );
		
		if(IndexLista(aux, lista, tamanhoLista) == -1)
			lista[cont++] = aux;
	}

	for(int i = 0; i < tamanhoLista; ++i)
	{
		fprintf(out, "%i\n", lista[i] );
	}

	free(lista);

	fclose(out);
}

void ListaSequencial(int limiteInferior, int tamanhoLista, char* nomeArq)
{
	FILE *out = stdout;
	
	out = fopen(nomeArq, "w");
	
	for(int i = 1; i <= tamanhoLista; ++i)
		fprintf(out, "%i\n", i );
	
	fclose(out);
}
