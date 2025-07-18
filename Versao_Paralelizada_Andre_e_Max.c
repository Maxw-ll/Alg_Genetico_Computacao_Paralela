// Universidade Federal do Maranhão
// Computação Paralela
// Alunos: André Luiz Ribeiro de Araujo Lima e Mawxell Pires Silva
// Versão Corrigida

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#define TRUE 1
#define FALSE 0

/*Código Funções de Fitness*/
#define F1 1
#define F2 2
#define F3 3
#define F4 4
#define F5 5

/*Quantidade Registradores Locais*/
#define NUM_REGS 4
/*Quantidade de Instruções que compõem o Cromossomo*/
#define NUM_INSTRUCOES 4
/*Cada Instrução tem 8 Bits -> 4 Bits para Representar o Código da Instrução e 4 Bits para Operandos (2 para cada Registrador) */
#define BITS_INSTRUCAO 8
/*Quantidade de Bits de um Registrador. Efetivamente o registrador é apenas um valor inteiro condicionado a valores no intervalor de 0 até 2^(n) - 1.*/
#define BITS_REG 8
/*O Cromossomo é composto apenas pelas instruções, que modificarão seus Registradores Locais*/
#define CHROMOSOME_LENGTH (NUM_INSTRUCOES * BITS_INSTRUCAO)
/*Tamanho da População*/
#define TAMPOP 1000
/*Quantidade Máxima de Gerações*/
#define MAXGER 100000 // Reduzido para testes
/*Taxa de Cruzamento*/
#define TAXCRUZ 0.8
/*Taxa de Mutação*/
#define TAXMUT 0.10
/*Elitismo*/
#define ELITISMO 1
/*Para Cromossomos em que ocorra Loop, ainda tentamos repetir até no Máximo 2 Vezes a Quantidade Original de Instruções em um Cromossomo*/
#define MAX_STEPS (NUM_INSTRUCOES * 2)
/*Fitness para Cromossomo que Ocorreu Loop*/
#define PENALIDADE_LOOP 999999

/*Tags para comunicação MPI*/
#define TAG_TRABALHO 1
#define TAG_RESULTADO 2
#define TAG_TERMINAR 3

/*Estrutura do Indivíduo*/
typedef struct
{
    unsigned char genes[CHROMOSOME_LENGTH];
    int regs[NUM_REGS];
    int fitness;
    int age_regs[NUM_REGS];
} Individuo;

/*Estrutura para enviar trabalho aos escravos*/
typedef struct
{
    unsigned char genes[CHROMOSOME_LENGTH];
    int regs[NUM_REGS];
    int tipo_fitness;
    int individuo_id;
    int age_regs[NUM_REGS];
} TrabalhoEscravo;

/*Estrutura para receber resultado dos escravos*/
typedef struct
{
    int fitness;
    int regs[NUM_REGS];
    int individuo_id;
    int age_regs[NUM_REGS];
} ResultadoEscravo;

/*População*/
Individuo pop[TAMPOP], nova_pop[TAMPOP];

/*Gera um Indivíduo Aleatório*/
void gera_individuo(Individuo *ind)
{
    for (int i = 0; i < CHROMOSOME_LENGTH; i++)
    {
        ind->genes[i] = rand() % 2;
    }

    for (int i = 0; i < NUM_REGS; i++)
    {
        ind->regs[i] = rand() % (int)pow(2, BITS_REG);
    }
    ind->fitness = 0;
}

/*Mutação*/
void mutacao(Individuo *ind)
{
    if (((float)rand() / RAND_MAX) < TAXMUT)
    {
        int mutacao = rand() % 3;

        switch (mutacao)
        {
        case 0:
            // Inverte alguns bits aleatórios (não todos)
            for (int i = 0; i < CHROMOSOME_LENGTH; i++)
            {
                if (((float)rand() / RAND_MAX) < 0.1) // 10% chance por bit
                {
                    ind->genes[i] ^= 1;
                }
            }
            break;
        case 1:
            // Sorteia uma Instrução e reescreve
            {
                int instr_index = rand() % NUM_INSTRUCOES;
                int base = instr_index * BITS_INSTRUCAO;
                for (int j = 0; j < BITS_INSTRUCAO; j++)
                {
                    ind->genes[base + j] = rand() % 2;
                }
            }
            break;
        case 2:
            // Troca duas instruções
            {
                int i1 = rand() % NUM_INSTRUCOES;
                int i2 = rand() % NUM_INSTRUCOES;
                int base_i1 = i1 * BITS_INSTRUCAO;
                int base_i2 = i2 * BITS_INSTRUCAO;
                for (int j = 0; j < BITS_INSTRUCAO; j++)
                {
                    unsigned char temp = ind->genes[base_i1 + j];
                    ind->genes[base_i1 + j] = ind->genes[base_i2 + j];
                    ind->genes[base_i2 + j] = temp;
                }
            }
            break;
        }
    }
}

/*Cruzamento*/
void cruzamento(Individuo pai1, Individuo pai2, Individuo *filho1, Individuo *filho2)
{
    int ponto = rand() % CHROMOSOME_LENGTH;
    for (int i = 0; i < CHROMOSOME_LENGTH; i++)
    {
        if (i < ponto)
        {
            filho1->genes[i] = pai1.genes[i];
            filho2->genes[i] = pai2.genes[i];
        }
        else
        {
            filho1->genes[i] = pai2.genes[i];
            filho2->genes[i] = pai1.genes[i];
        }
    }
}

/*Calcula o Fitness*/
int calcula_fitness(int *regs, int tipo)
{
    int resultado = 0;
    switch (tipo)
    {
    case F1: /* A + B = C + D*/
        resultado = abs((regs[0] + regs[1]) - (regs[2] + regs[3]));
        break;
    case F2: /* A + B > C - D*/
        if ((regs[0] + regs[1]) > (regs[2] - regs[3]))
        {
            resultado = 0;
        }
        else
        {
            resultado = abs((regs[0] + regs[1]) - (regs[2] - regs[3]));
        }
        break;
    case F3: /* A % B == 0 -> B divide A*/
        if (regs[1] != 0 && regs[0] % regs[1] == 0)
        {
            resultado = 0;
        }
        else
        {
            if (regs[1] == 0)
                resultado = (int)pow(2, BITS_REG) - 1;
            else
                resultado = regs[0] % regs[1];
        }
        break;
    case F4: /* A=B+1 && B=C+1 && C=D+1 */
    {
        int diff1 = abs(regs[0] - (regs[1] + 1));
        int diff2 = abs(regs[1] - (regs[2] + 1));
        int diff3 = abs(regs[2] - (regs[3] + 1));
        resultado = diff1 + diff2 + diff3;
    }
    break;
    case F5: /*CH(A)+CH(B)+CH(C)+CH(D)+CH(A) = "OTIMO"*/
    {
        char alvo[] = {'O', 'T', 'I', 'M', 'O'};
        int match = 0;
        int valores[5] = {regs[0], regs[1], regs[2], regs[3], regs[0]};
        for (int i = 0; i < 5; i++)
        {
            if ((char)valores[i] == alvo[i])
            {
                match++;
            }
        }
        resultado = 5 - match;
    }
    break;
    }
    return resultado;
}

/*Normalizar registradores*/
void normalized_regs(int *regs)
{
    for (int i = 0; i < NUM_REGS; i++)
    {
        if (regs[i] < 0)
        {
            regs[i] = 0;
        }
        else if (regs[i] >= (int)pow(2, BITS_REG))
        {
            regs[i] = (int)pow(2, BITS_REG) - 1;
        }
    }
}

/*Executa uma Instrução*/
int executa_instrucao(unsigned char *instrucao, int *regs, int pc)
{

    int opcode = (instrucao[0] << 3) | (instrucao[1] << 2) | (instrucao[2] << 1) | instrucao[3];
    int r1 = (instrucao[4] << 1) | instrucao[5];
    int r2 = (instrucao[6] << 1) | instrucao[7];

    switch (opcode)
    {
    case 0: // ADD
        regs[r1] += regs[r2];
        break;
    case 1: // SUB
        regs[r1] -= regs[r2];
        break;
    case 2: // MULT
        regs[r1] *= regs[r2];
        break;
    case 3: // DIV
        if (regs[r2] != 0)
        {
            regs[r1] /= regs[r2];
        }
        break;
    case 4: // MOD
        if (regs[r2] != 0)
        {
            regs[r1] %= regs[r2];
        }
        break;
    case 5: // INC
        regs[r1]++;
        break;
    case 6: // DEC
        regs[r1]--;
        break;
    case 7: // MOV
        regs[r1] = regs[r2];
        break;
    case 8: // NOT
        regs[r1] = ~regs[r1];
        break;
    case 9: // GT
        regs[r1] = (regs[r1] > regs[r2]) ? TRUE : FALSE;
        break;
    case 10: // LT
        regs[r1] = (regs[r1] < regs[r2]) ? TRUE : FALSE;
        break;
    case 11: // EQ
        regs[r1] = (regs[r1] == regs[r2]) ? TRUE : FALSE;
        break;
    case 12: // JMP
        return regs[r1] % NUM_INSTRUCOES;
    case 13: // JZ
        return (regs[r1] == 0) ? regs[r2] % NUM_INSTRUCOES : pc + 1;
    case 14: // JNZ
        return (regs[r1] != 0) ? regs[r2] % NUM_INSTRUCOES : pc + 1;
    case 15: // NOP
        break;
    }
    normalized_regs(regs);
    return pc + 1;
}

/*Avaliação paralela*/
int avaliar_paralelo(unsigned char *genes, int *regs, int *age_regs, int tipo_fitness)
{
    int pc = 0, steps = 0;
    int regs_temp[NUM_REGS];

    for (int i = 0; i < NUM_REGS; i++)
    {
        regs_temp[i] = regs[i];
    }

    while (steps < MAX_STEPS && pc < NUM_INSTRUCOES)
    {
        pc = executa_instrucao(&genes[pc * BITS_INSTRUCAO], regs_temp, pc);
        steps++;
    }

    if (steps >= MAX_STEPS)
    {
        return PENALIDADE_LOOP;
    }

    // Copiar registradores finais de volta
    for (int i = 0; i < NUM_REGS; i++)
    {
        age_regs[i] = regs[i];
    }

    for (int i = 0; i < NUM_REGS; i++)
    {
        regs[i] = regs_temp[i];
    }

    return calcula_fitness(regs_temp, tipo_fitness);
}

void salvar_melhor_individuo(int idx)
{
    const char *ops[] = {
        "ADD", "SUB", "MUL", "DIV", "MOD", "INC", "DEC", "MOV",
        "NOT", "GT", "LT", "EQ", "JMP", "JZ", "JNZ", "NOP"};

    FILE *f = fopen("Melhor_Individuo_MPI.txt", "w");

    Individuo *ind = &pop[idx];

    fprintf(f, "Individuo %d \n", idx);

    fprintf(f, "Cromossomo: ");

    for (int i = 0; i < CHROMOSOME_LENGTH; i++)
    {
        fprintf(f, "%d", ind->genes[i]);
    }
    fprintf(f, "\n");

    fprintf(f, "Registradores -> Valores Iniciais: ");
    for (int i = 0; i < NUM_REGS; i++)
    {
        fprintf(f, "R%d=%d ", i, ind->age_regs[i]);
    }

    fprintf(f, "\n");

    fprintf(f, "Instrucoes:\n");

    int pc = 0, steps = 0;

    while (steps < MAX_STEPS && pc < NUM_INSTRUCOES)
    {
        int pos = pc * BITS_INSTRUCAO;
        int opcode = (ind->genes[pos] << 3) | (ind->genes[pos + 1] << 2) | (ind->genes[pos + 2] << 1) | ind->genes[pos + 3];
        int r1 = (ind->genes[pos + 4] << 1) | ind->genes[pos + 5];
        int r2 = (ind->genes[pos + 6] << 1) | ind->genes[pos + 7];

        fprintf(f, "%s R%d, R%d\n", ops[opcode], r1, r2);
        steps++;
        pc++;
    }

    // zerar_pcs(ind);

    fprintf(f, "Registradores -> Valores Finais: ");
    for (int i = 0; i < NUM_REGS; i++)
    {
        fprintf(f, "R%d=%d ", i, ind->regs[i]);
    }
    fprintf(f, "\n");
    fprintf(f, "Fitness Final: %d\n", ind->fitness);
    fprintf(f, "\n--------------------------\n");

    fclose(f);
}

/*Avaliação da população paralela - CORRIGIDA*/
void avaliar_populacao_paralela(int tipo_fitness, int rank, int size)
{
    if (rank == 0) // Processo Mestre
    {
        TrabalhoEscravo trabalho;
        ResultadoEscravo resultado;
        int enviados = 0;
        int recebidos = 0;

        // Enviar trabalho inicial
        for (int escravo = 1; escravo < size && enviados < TAMPOP; escravo++)
        {
            trabalho.tipo_fitness = tipo_fitness;
            trabalho.individuo_id = enviados;
            memcpy(trabalho.genes, pop[enviados].genes, CHROMOSOME_LENGTH);
            memcpy(trabalho.regs, pop[enviados].regs, NUM_REGS * sizeof(int));
            memcpy(trabalho.age_regs, pop[enviados].age_regs, NUM_REGS * sizeof(int));

            MPI_Send(&trabalho, sizeof(TrabalhoEscravo), MPI_BYTE, escravo, TAG_TRABALHO, MPI_COMM_WORLD);
            enviados++;
        }

        // Receber resultados e enviar novos trabalhos
        while (recebidos < TAMPOP)
        {
            MPI_Status status;
            MPI_Recv(&resultado, sizeof(ResultadoEscravo), MPI_BYTE, MPI_ANY_SOURCE, TAG_RESULTADO, MPI_COMM_WORLD, &status);

            // Aplicar resultado
            pop[resultado.individuo_id].fitness = resultado.fitness;
            memcpy(pop[resultado.individuo_id].regs, resultado.regs, NUM_REGS * sizeof(int));
            memcpy(pop[resultado.individuo_id].age_regs, resultado.age_regs, NUM_REGS * sizeof(int));
            recebidos++;

            // Enviar novo trabalho se houver
            if (enviados < TAMPOP)
            {
                trabalho.tipo_fitness = tipo_fitness;
                trabalho.individuo_id = enviados;
                memcpy(trabalho.genes, pop[enviados].genes, CHROMOSOME_LENGTH);
                memcpy(trabalho.regs, pop[enviados].regs, NUM_REGS * sizeof(int));
                memcpy(trabalho.age_regs, pop[enviados].age_regs, NUM_REGS * sizeof(int));

                MPI_Send(&trabalho, sizeof(TrabalhoEscravo), MPI_BYTE, status.MPI_SOURCE, TAG_TRABALHO, MPI_COMM_WORLD);
                enviados++;
            }
        }
    }
    else // Processos Escravos
    {
        TrabalhoEscravo trabalho;
        ResultadoEscravo resultado;
        MPI_Status status;

        while (1)
        {
            MPI_Recv(&trabalho, sizeof(TrabalhoEscravo), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (status.MPI_TAG == TAG_TERMINAR)
            {
                break;
            }

            // Processar trabalho
            int regs_temp[NUM_REGS];
            int age_regs_temps[NUM_REGS];
            memcpy(regs_temp, trabalho.regs, NUM_REGS * sizeof(int));
            memcpy(age_regs_temps, trabalho.age_regs, NUM_REGS * sizeof(int));

            resultado.fitness = avaliar_paralelo(trabalho.genes, regs_temp, age_regs_temps,trabalho.tipo_fitness);
            resultado.individuo_id = trabalho.individuo_id;
            memcpy(resultado.regs, regs_temp, NUM_REGS * sizeof(int));
            memcpy(resultado.age_regs, age_regs_temps, NUM_REGS * sizeof(int));

            MPI_Send(&resultado, sizeof(ResultadoEscravo), MPI_BYTE, 0, TAG_RESULTADO, MPI_COMM_WORLD);
        }
    }
}

/*Torneio*/
int torneio()
{
    int a = rand() % TAMPOP;
    int b = rand() % TAMPOP;
    return (pop[a].fitness < pop[b].fitness) ? a : b;
}

/*Copiar indivíduo*/
void copiar_individuo(Individuo *src, Individuo *dst)
{
    memcpy(dst->genes, src->genes, CHROMOSOME_LENGTH);
    memcpy(dst->regs, src->regs, NUM_REGS * sizeof(int));
    dst->fitness = src->fitness;
}

/*Gerar novas entradas - MODIFICADO para ser menos agressivo*/
void gerar_novas_entradas()
{
    for (int i = 0; i < TAMPOP; i++)
    {
        // Só muda os registradores com 30% de chance
        if (((float)rand() / RAND_MAX) < 0.3)
        {
            for (int j = 0; j < NUM_REGS; j++)
            {
                if (rand() % 2 == 0)
                {
                    pop[i].regs[j] += (rand() % 5) + 1;
                }
                else
                {
                    pop[i].regs[j] -= (rand() % 5) + 1;
                }
            }
            normalized_regs(pop[i].regs);
        }
    }
}

/*Nova geração*/
void nova_geracao()
{
    int i = 0;
    int melhor = 0;

    for (int j = 1; j < TAMPOP; j++)
    {
        if (pop[j].fitness < pop[melhor].fitness)
            melhor = j;
    }

    if (ELITISMO)
    {
        copiar_individuo(&pop[melhor], &nova_pop[i++]);
    }

    while (i < TAMPOP)
    {
        int p1 = torneio();
        int p2 = torneio();
        Individuo f1, f2;

        if (((float)rand() / RAND_MAX) < TAXCRUZ)
        {
            cruzamento(pop[p1], pop[p2], &f1, &f2);
        }
        else
        {
            copiar_individuo(&pop[p1], &f1);
            copiar_individuo(&pop[p2], &f2);
        }

        mutacao(&f1);
        mutacao(&f2);

        nova_pop[i] = f1;
        i++;
        if (i < TAMPOP)
        {
            nova_pop[i] = f2;
            i++;
        }
    }

    for (int i = 0; i < TAMPOP; i++)
    {
        copiar_individuo(&nova_pop[i], &pop[i]);
    }

    gerar_novas_entradas();
}

/*Finalizar escravos - CORRIGIDO*/
void finalizar_escravos(int size)
{
    int dummy = 0;
    for (int i = 1; i < size; i++)
    {
        MPI_Send(&dummy, sizeof(int), MPI_BYTE, i, TAG_TERMINAR, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[])
{
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL) + rank);

    int tipo_fitness = F5;

    if (rank == 0) // Processo Mestre
    {
        printf("Executando com %d processos (1 mestre + %d escravos)\n", size, size - 1);
        printf("Função de fitness: F%d\n", tipo_fitness);

        // Gerar população inicial
        for (int i = 0; i < TAMPOP; i++)
        {
            gera_individuo(&pop[i]);
        }

        // Avaliar população inicial
        avaliar_populacao_paralela(tipo_fitness, rank, size);

        FILE *f = fopen("Log_Geracoes_MPI.txt", "w");

        for (int g = 0; g < MAXGER; g++)
        {
            int melhor = 0;
            int pior = 0;
            double media = 0.0;

            for (int i = 0; i < TAMPOP; i++)
            {
                if (pop[i].fitness < pop[melhor].fitness)
                    melhor = i;
                if (pop[i].fitness > pop[pior].fitness)
                    pior = i;
                media += pop[i].fitness;
            }
            media /= TAMPOP;

            printf("Ger %d | Melhor: %d\n", g, pop[melhor].fitness);
            fprintf(f, "Ger %d | Melhor: %d\n", g, pop[melhor].fitness);

            // Critério de parada
            if (pop[melhor].fitness == 0)
            {
                printf("Solução ótima encontrada na geração %d!\n", g);
                salvar_melhor_individuo(melhor);

                break;
            }

            nova_geracao();
            avaliar_populacao_paralela(tipo_fitness, rank, size);
        }

        fclose(f);

        // Finalizar processos escravos
        finalizar_escravos(size);

        // Exibir melhor indivíduo
        int melhor_final = 0;
        for (int i = 1; i < TAMPOP; i++)
        {
            if (pop[i].fitness < pop[melhor_final].fitness)
                melhor_final = i;
        }

        printf("\nMelhor indivíduo final:\n");
        printf("Fitness: %d\n", pop[melhor_final].fitness);
        printf("Registradores: ");
        for (int i = 0; i < NUM_REGS; i++)
        {
            printf("R%d=%d ", i, pop[melhor_final].regs[i]);
        }
        printf("\n");
    }
    else // Processos Escravos
    {
        // Escravos processam todas as gerações
        for (int g = 0; g <= MAXGER; g++)
        {
            avaliar_populacao_paralela(tipo_fitness, rank, size);
        }
    }

    MPI_Finalize();
    return 0;
}