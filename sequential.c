#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 

// Estados possíveis
#define HEALTH 1
#define INFECTED -1
#define DEAD -2
#define EMPTY 0

int hasContaminatingNeighbor(int **matrix, int i, int j, int N, int M);

int main(int argc, char **argv)
{
    if(argc < 2)
    {
        perror("You must provide a .txt file.");
        return 1;
    }

    FILE *input = fopen(argv[1], "r");
    if(!input)
    {
        perror("Can't open file");
        return 1;
    }

    int N, M;
    fscanf(input, "%d %d", &N, &M);

    int **matrix = malloc(N * sizeof(int *));
    int **next = malloc(N * sizeof(int *));
    int **dead_age = malloc(N * sizeof(int *));
    for(int i = 0; i < N; i++)
    {
        matrix[i] = malloc(M * sizeof(int));
        next[i] = malloc(M * sizeof(int));
        dead_age[i] = calloc(M, sizeof(int)); // Inicializa já com 0. 
    }

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            fscanf(input, "%d", &matrix[i][j]);
        }
    }
    fclose(input);

    // Caso haja mortos desde a primeira matriz. 
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            if(matrix[i][j] == DEAD)
                dead_age[i][j] = 1;
        }
    }



    srand(time(NULL));
    int max_iter = N * M; 
    int iter = 0;
   
    // ------------ INÍCIO DA CAPTURA DE TEMPO ------------------------------
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    while(iter < max_iter)
    {
        int living = 0, dead = 0;
        int flag = 0; // Flag para quando uma pessoa mudou de estado (saudável, morto ou infectado).

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < M; j++)
            {
                int current = matrix[i][j];

                if(current == HEALTH) // Se for saudável, procura por um vizinho contaminante.
                {
                    if(hasContaminatingNeighbor(matrix, i, j, N, M)) // Possui vizinho contaminante, então se torna infectado. 
                    {
                        next[i][j] = INFECTED; 
                        flag = 1;
                    }
                    else // Não possui vizinho contaminante, segue saudável. 
                        next[i][j] = HEALTH;
                }
                else if(current == INFECTED)
                {
                    int probability = rand() %  10000; 

                    if(probability <= 999) next[i][j] = HEALTH; // curado
                    else if(probability <= 3999) next[i][j] = INFECTED; // continua infectado
                    else 
                    {    
                        next[i][j] = DEAD; // morreu
                        dead_age[i][j] = 1;
                    }
                    flag = 1;
                }
                else if(current == DEAD)
                {
                    if(dead_age[i][j] == 1)
                    {
                        next[i][j] = DEAD; // Continua morto por mais uma iteração.
                        dead_age[i][j] = 2; // Aumenta a "idade" do morto. 
                    }
                    else if(dead_age[i][j] == 2)
                    {
                        next[i][j] = EMPTY; // Morto sumiu depois de duas iterações. 
                        dead_age[i][j] = 0; // Zera a idade dele. 
                        flag = 1;
                    }
                }
                else
                    next[i][j] = EMPTY;
                
                if(next[i][j] == HEALTH || next[i][j] == INFECTED)
                    living++;
                else if(next[i][j] == DEAD)
                    dead++;
            }   
        }

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < M; j++)
            {
                matrix[i][j] = next[i][j];
            }
        }

        iter++;

        // Condição de parada antecipada: ninguém foi alterado ou não há mais vivos. 
        if(!flag || living == 0)
            break;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    // ------------------------- FIM DA CAPTURA DE TEMPO -----------------------------------------
    // converte para milissegundos
    double elapsed_ms = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1e6;
    printf("Tempo CPU: %.3f ms\n", elapsed_ms);

    int totalDeaths = 0, totalAlive = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            if(matrix[i][j] == DEAD) totalDeaths++;
            else if(matrix[i][j] == HEALTH || matrix[i][j] == INFECTED) totalAlive++;
        }
    }

    FILE *output = fopen("output.txt", "w");
    fprintf(output, "Mortos: %d \n Sobreviventes: %d \n", totalDeaths, totalAlive);
    fclose(output);

    for(int i = 0; i < N; i++)
    {
            free(matrix[i]);
            free(next[i]);
    }

    free(matrix);
    free(next);

    return 0;
}   

int hasContaminatingNeighbor(int **matrix, int i, int j, int N, int M)
{
    /*
        Todas as direções possíveis: 
        {-1, 0} = cima;
        {1, 0} = baixo;
        {0, -1} = esquerda;
        {0, 1} = direita;
    */
    int directions[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    for(int d = 0; d < 4; d++)
    {
        // Calcula a posição de cada um dos 4 possíveis vizinhos. 
        int ni = i + directions[d][0];
        int nj = j + directions[d][1];

        // Verifica se o possível vizinho realmente existe. 
        if(ni >= 0 && ni < N && nj >= 0 && nj < M)
        {
            if(matrix[ni][nj] == INFECTED || matrix[ni][nj] == DEAD)
                return 1;
        }
    }

    return 0;
}