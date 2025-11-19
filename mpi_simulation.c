#define _POSIX_C_SOURCE 199309L
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Estados possíveis
#define HEALTH 1
#define INFECTED -1
#define DEAD -2
#define EMPTY 0

// Retorna 1 se a célula (local_i, j) (considerando ghost rows) possui vizinho contaminante
int hasContaminatingNeighbor_local(int **local_matrix, int local_i, int j, int local_rows_plus_2, int M)
{
    // local_i está em [0 .. local_rows_plus_2-1], onde 0 e local_rows_plus_2-1 são ghost rows.
    // Vizinhança: cima, baixo, esquerda, direita (apenas horizontal/vertical).
    int directions[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for(int d = 0; d < 4; d++)
    {
        int ni = local_i + directions[d][0];
        int nj = j + directions[d][1];

        if(ni >= 0 && ni < local_rows_plus_2 && nj >= 0 && nj < M)
        {
            if(local_matrix[ni][nj] == INFECTED || local_matrix[ni][nj] == DEAD)
                return 1;
        }
    }
    return 0;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(argc < 2)
    {
        if(rank == 0) perror("You must provide a .txt file.");
        MPI_Finalize();
        return 1;
    }

    int N = 0, M = 0;
    int *full_matrix_flat = NULL; // apenas em rank 0: N*M ints lidos em ordem row-major

    if(rank == 0)
    {
        FILE *input = fopen(argv[1], "r");
        if(!input)
        {
            perror("Can't open file");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
        if(fscanf(input, "%d %d", &N, &M) != 2)
        {
            perror("Error reading N and M");
            fclose(input);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        // aloca e lê a matriz inteira em flat array
        full_matrix_flat = malloc(N * M * sizeof(int));
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < M; j++)
            {
                int v;
                if(fscanf(input, "%d", &v) != 1) v = 0;
                full_matrix_flat[i * M + j] = v;
            }
        }
        fclose(input);
    }

    // Broadcast de N e M para todos os processos
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(N <= 0 || M <= 0)
    {
        if(rank == 0) fprintf(stderr, "N and M must be positive.\n");
        if(full_matrix_flat) free(full_matrix_flat);
        MPI_Finalize();
        return 1;
    }

    // Calcula quantas linhas cada processo recebe (balanceamento com resto)
    int base = N / nprocs;
    int rem = N % nprocs;
    int *rows_per_proc = malloc(nprocs * sizeof(int));
    int *displs_rows = malloc(nprocs * sizeof(int)); // displacement em linhas
    for(int p = 0; p < nprocs; p++)
    {
        rows_per_proc[p] = base + (p < rem ? 1 : 0);
    }
    displs_rows[0] = 0;
    for(int p = 1; p < nprocs; p++)
        displs_rows[p] = displs_rows[p-1] + rows_per_proc[p-1];

    int local_rows = rows_per_proc[rank]; // número de linhas "reais" deste processo
    int recvcount = local_rows * M;

    // Buffer para receber a parte flat (somente linhas reais, sem ghosts)
    int *local_flat = malloc((local_rows * M) * sizeof(int));
    // Prepare sendcounts and displs (em elementos, para MPI_Scatterv)
    int *sendcounts = NULL;
    int *senddispls = NULL;
    if(rank == 0)
    {
        sendcounts = malloc(nprocs * sizeof(int));
        senddispls = malloc(nprocs * sizeof(int));
        for(int p = 0; p < nprocs; p++)
        {
            sendcounts[p] = rows_per_proc[p] * M;
            senddispls[p] = displs_rows[p] * M;
        }
    }

    // Scatterv para distribuir as linhas em row-major contiguous
    MPI_Scatterv(full_matrix_flat, sendcounts, senddispls, MPI_INT,
                 local_flat, recvcount, MPI_INT,
                 0, MPI_COMM_WORLD);

    // libere full_matrix_flat no rank0, pois já não é mais necessário
    if(rank == 0 && full_matrix_flat)
    {
        free(full_matrix_flat);
        full_matrix_flat = NULL;
    }

    // Criar matrizes locais com ghost rows:
    // estrutura: linhas [0 .. local_rows+1], onde 0 = ghost topo, local_rows+1 = ghost base.
    int local_rows_plus_2 = local_rows + 2;
    int **local_matrix = malloc(local_rows_plus_2 * sizeof(int *));
    int **local_next = malloc(local_rows_plus_2 * sizeof(int *));
    int **local_dead_age = malloc(local_rows_plus_2 * sizeof(int *));
    for(int i = 0; i < local_rows_plus_2; i++)
    {
        local_matrix[i] = malloc(M * sizeof(int));
        local_next[i] = malloc(M * sizeof(int));
        local_dead_age[i] = calloc(M, sizeof(int)); // inicializa com zeros
    }

    // Preencher linhas reais (1..local_rows) a partir de local_flat
    for(int i = 0; i < local_rows; i++)
    {
        for(int j = 0; j < M; j++)
        {
            int val = local_flat[i * M + j];
            local_matrix[i+1][j] = val;
            // inicializa dead_age se já houver mortos na entrada
            if(val == DEAD)
                local_dead_age[i+1][j] = 1;
            else
                local_dead_age[i+1][j] = 0;
            local_next[i+1][j] = EMPTY;
        }
    }

    // Inicializa ghost rows (0 e local_rows+1) como EMPTY
    for(int j = 0; j < M; j++)
    {
        local_matrix[0][j] = EMPTY;
        local_matrix[local_rows+1][j] = EMPTY;
        local_next[0][j] = EMPTY;
        local_next[local_rows+1][j] = EMPTY;
        local_dead_age[0][j] = 0;
        local_dead_age[local_rows+1][j] = 0;
    }

    free(local_flat);
    if(rank == 0)
    {
        free(sendcounts);
        free(senddispls);
    }

    // Prepara comunicação entre vizinhos
    int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int down = (rank == nprocs - 1) ? MPI_PROC_NULL : rank + 1;

    // parâmetros da simulação
    int max_iter = N * M;
    int iter = 0;

    // semente rand por processo
    srand(time(NULL) + rank * 7919);

    // Barrier e medição do tempo
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    while(iter < max_iter)
    {
        // Troca de bordas: enviamos nossa primeira real linha (índice 1) para up, e recebemos a borda superior em ghost 0.
        // Enviamos nossa última real linha (índice local_rows) para down, e recebemos a borda inferior em ghost local_rows+1.
        // Usamos MPI_Sendrecv para evitar deadlocks.
        // Top ghost (row 0) <- last row of up
        MPI_Sendrecv(local_matrix[1], M, MPI_INT, up, 0,
                     local_matrix[0], M, MPI_INT, up, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Bottom ghost (row local_rows+1) <- first row of down
        MPI_Sendrecv(local_matrix[local_rows], M, MPI_INT, down, 1,
                     local_matrix[local_rows+1], M, MPI_INT, down, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Observação: se up/down == MPI_PROC_NULL, MPI_Sendrecv é noop apropriado.

        // --- Computação local: aplica as regras para linhas reais [1 .. local_rows] ---
        int flag_local = 0;   // se alguma célula mudou
        int living_local = 0; // contagem de vivos (health ou infected) na próxima matriz
        int dead_local = 0;   // contagem de mortos na próxima matriz

        for(int i = 1; i <= local_rows; i++)
        {
            for(int j = 0; j < M; j++)
            {
                int current = local_matrix[i][j];

                if(current == HEALTH)
                {
                    if(hasContaminatingNeighbor_local(local_matrix, i, j, local_rows_plus_2, M))
                    {
                        local_next[i][j] = INFECTED;
                        flag_local = 1;
                    }
                    else
                        local_next[i][j] = HEALTH;
                }
                else if(current == INFECTED)
                {
                    int probability = rand() % 10000;
                    if(probability <= 999)
                    {
                        local_next[i][j] = HEALTH; // curado
                    }
                    else if(probability <= 3999)
                    {
                        local_next[i][j] = INFECTED; // continua infectado
                    }
                    else
                    {
                        local_next[i][j] = DEAD; // morreu
                        // define dead_age na próxima iteração: setamos agora para 1 para indicar que recém morreu
                        // porém dead_age pertence ao "estado atual" e precisa ser usado nas iterações seguintes.
                        // vamos setar dead age temporariamente; será copiado ao atualizar.
                        local_dead_age[i][j] = 1;
                    }
                    flag_local = 1;
                }
                else if(current == DEAD)
                {
                    if(local_dead_age[i][j] == 1)
                    {
                        // Continua morto por mais uma iteração.
                        local_next[i][j] = DEAD;
                        // incrementa idade do morto
                        // mas incrementaremos local_dead_age após copiar next->matrix
                        // setamos para 2 agora
                        local_dead_age[i][j] = 2;
                    }
                    else if(local_dead_age[i][j] == 2)
                    {
                        // Morto some depois de duas iterações.
                        local_next[i][j] = EMPTY;
                        local_dead_age[i][j] = 0;
                        flag_local = 1;
                    }
                    else
                    {
                        // Caso estranho (DEAD mas dead_age==0) — mantemos como DEAD e setamos idade para 1
                        local_next[i][j] = DEAD;
                        local_dead_age[i][j] = 1;
                        flag_local = 1;
                    }
                }
                else // EMPTY
                {
                    local_next[i][j] = EMPTY;
                }

                if(local_next[i][j] == HEALTH || local_next[i][j] == INFECTED)
                    living_local++;
                else if(local_next[i][j] == DEAD)
                    dead_local++;
            }
        }

        // Precisamos também atualizar dead_age corretamente:
        // As regras do seu sequencial eram:
        // - quando alguém morre, dead_age = 1
        // - se dead_age == 1 -> na próxima iteração permanece DEAD e passa a dead_age 2
        // - se dead_age == 2 -> some (EMPTY) e dead_age=0
        // Em nosso laço acima já atualizamos local_dead_age nos casos correspondentes.
        // Agora copiamos local_next -> local_matrix para linhas reais e também tratamos ghost rows.

        for(int i = 1; i <= local_rows; i++)
        {
            for(int j = 0; j < M; j++)
            {
                // É importante resetar dead_age para células que não estão mortas.
                if(local_next[i][j] == DEAD)
                {
                    // Se já estava DEAD e pizza local_dead_age foi atualizado (1->2), mantemos o valor.
                    // Se foi recém-morto (pela transição INFECTED->DEAD) e definimos local_dead_age[i][j]=1 acima.
                    // Caso local_dead_age esteja 0 por algum motivo, setamos 1.
                    if(local_dead_age[i][j] == 0)
                        local_dead_age[i][j] = 1;
                }
                else
                {
                    // não morto -> zera dead_age
                    local_dead_age[i][j] = 0;
                }

                local_matrix[i][j] = local_next[i][j];
            }
        }

        // Após atualizar, precisamos checar se houve mudança global: usamos MPI_Allreduce (OR lógico).
        int flag_global = 0;
        MPI_Allreduce(&flag_local, &flag_global, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

        // Também precisamos saber se há vivos no mundo: reduzir living_local por soma.
        int living_global = 0;
        MPI_Allreduce(&living_local, &living_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        iter++;

        if(!flag_global || living_global == 0)
            break;
    } // fim do while

    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    double elapsed_ms = (end - start) * 1000.0;

    // Contagens finais: cada processo conta mortos e vivos em suas linhas reais
    int local_dead_final = 0;
    int local_alive_final = 0;
    for(int i = 1; i <= local_rows; i++)
    {
        for(int j = 0; j < M; j++)
        {
            if(local_matrix[i][j] == DEAD) local_dead_final++;
            else if(local_matrix[i][j] == HEALTH || local_matrix[i][j] == INFECTED) local_alive_final++;
        }
    }

    int global_dead = 0, global_alive = 0;
    MPI_Reduce(&local_dead_final, &global_dead, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_alive_final, &global_alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // O processo 0 escreve output e imprime tempo
    if(rank == 0)
    {
        FILE *output = fopen("output.txt", "w");
        if(output)
        {
            fprintf(output, "Mortos: %d \nSobreviventes: %d \n", global_dead, global_alive);
            fclose(output);
        }
        printf("Tempo (ms): %.3f\n", elapsed_ms);
    }

    // libera memória
    for(int i = 0; i < local_rows_plus_2; i++)
    {
        free(local_matrix[i]);
        free(local_next[i]);
        free(local_dead_age[i]);
    }
    free(local_matrix);
    free(local_next);
    free(local_dead_age);

    free(rows_per_proc);
    free(displs_rows);

    MPI_Finalize();
    return 0;
}
