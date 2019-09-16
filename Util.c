void Limites() {

   printf("MIN SIGNED CHAR = %d\n", SCHAR_MIN);
   printf("MAX SIGNED CHAR = %d\n", SCHAR_MAX);
   printf("MAX UNSIGNED CHAR = %d\n", UCHAR_MAX);

   printf("MIN SHORT INT = %d\n", SHRT_MIN);
   printf("MAX SHORT INT = %d\n", SHRT_MAX); 

   printf("MIN INT = %d\n", INT_MIN);
   printf("MAX INT = %d\n", INT_MAX);

   printf("MIN CHAR = %d\n", CHAR_MIN);
   printf("MAX CHAR = %d\n", CHAR_MAX);

   printf("MIN LONG = %ld\n", LONG_MIN);
   printf("MAX LONG = %ld\n", LONG_MAX);
}

//============================================================================================================

int IndexLista(int item, int *lista, int sizeLista)
{
	int ind = -1;

	for(int i = 0; i < sizeLista; ++i)
	{
		if(lista[i] == item)
		{
			ind = i;
			break;
		}
	}

	return ind;
}

//============================================================================================================

void PrintLista(int *lista, int sizeLista)
{
	for(int i = 0; i < sizeLista; ++i)
	{
		printf("%i ", lista[i]);
	}

	printf("\n");
}

//============================================================================================================

void PrintMatriz(int **matriz, int sizeMatriz, int sizeLista)
{
	printf("Matriz [%i][%i]\n", sizeMatriz, sizeLista);

	for(int i = 0; i < sizeMatriz; ++i)
	{
		for(int j = 0; j < sizeLista; ++j)
			printf("%i ", matriz[i][j]);

		printf("\n");
	}

	printf("\n");
}

