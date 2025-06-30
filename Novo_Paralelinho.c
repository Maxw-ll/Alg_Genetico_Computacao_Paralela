// Universidade Federal do Maranhão
// Computação Paralela
// Alunos: André Luiz Ribeiro de Araujo Lima e Mawxell Pires Silva

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

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
#define TAMPOP 5
/*Quantidade Máxima de Gerações*/
#define MAXGER 20
/*Taxa de Cruzamento*/
#define TAXCRUZ 0.8
/*Taxa de Mutação*/
#define TAXMUT 0.05
/*Elitismo*/
#define ELITISMO 1
/*Para Cromossomos em que ocorra Loop, ainda tentamos repetir até no Máximo 2 Vezes a Quantidade Original de Instruções em um Cromossomo*/
#define MAX_STEPS (NUM_INSTRUCOES * 2)
/*Fitness para Cromossomo que Ocorreu Loop*/
#define PENALIDADE_LOOP 999999

/*Estrutura do Indivíduo Composto pelo Cromossomo (A Parte que é evoluída) e seus Registradores, além de sua aptidão*/
typedef struct
{
    unsigned char genes[CHROMOSOME_LENGTH];
    int regs[NUM_REGS];
    int fitness;
} Individuo;

/*Inicia uma população de Indivíduos */
Individuo pop[TAMPOP], nova_pop[TAMPOP];

/*Gera um Indivíduo Aleatório -> Com Bits todos Aleatórios*/
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

/*Se um número aleatório entre 0 e 1 for menor que a taxa de mutação, então um certo indíviduo sofre mutação*/
void mutacao(Individuo *ind)
{
    if (((float)rand() / RAND_MAX) < TAXMUT)
    {
        int mutacao = rand() % 3;

        switch (mutacao)
        {
        case 0:
        {
            /*Inverte os Bits*/
            for (int i = 0; i < CHROMOSOME_LENGTH; i++)
            {
                ind->genes[i] ^= 1;
            }
        }
        break;
        case 1:
            /*Sorteia uma Instrução e reescreve de forma aleatória*/
            int instr_index = rand() % NUM_INSTRUCOES;
            int base = instr_index * BITS_INSTRUCAO;
            for (int j = 0; j < BITS_INSTRUCAO; j++)
            {
                ind->genes[base + j] = rand() % 2;
            }
            break;
        case 2:

            /*Troca uma Instrução com outra, escolhidas aleatoriamente. Pode ser sorteada a mesma*/
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
    }
}

/*Realiza o Cruzamento entre dois pais, misturando seus códigos genéticos entre um ponto definido aleatoriamente*/
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

/*Calcula o Fitness Dependendo da Função Desejada*/
int calcula_fitness(int *regs, int tipo)
{
    int resultado = 0;
    switch (tipo)
    {
    case 1: /* A + B = C + D*/
        resultado = abs((regs[0] + regs[1]) - (regs[2] + regs[3]));
        break;
    case 2: /* A + B > C - D*/
        if ((regs[0] + regs[1]) > regs[2] - regs[3])
        {
            resultado = TRUE;
        }
        else
        {
            resultado = FALSE;
        }
        break;
    case 3: /* A % B == 0 -> B divide A*/
        if (regs[1] != 0 && regs[0] % regs[1] == 0)
        {
            resultado = TRUE;
        }
        else
        {
            resultado = FALSE;
        }
        break;
    case 4: /* A=B+1 && B=C+1 && C=D+1 */
        if (regs[0] == regs[1] + 1 && regs[1] == regs[2] + 1 && regs[2] == regs[3] + 1)
        {
            resultado = TRUE;
        }
        else
        {
            resultado = FALSE;
        }
        break;
    case 5: /*CH(A)+CH(B)+CH(C)+CH(D)+CH(A) = "OTIMO"*/
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
        break;
    }
    }
    return resultado;
}

int executa_instrucao(unsigned char *instrucao, int *regs, int pc)
{
    int opcode = (instrucao[0] << 3) | (instrucao[1] << 2) | (instrucao[2] << 1) | instrucao[3];
    int r1 = (instrucao[4] << 1) | instrucao[5];
    int r2 = (instrucao[6] << 1) | instrucao[7];

    opcode %= 16;

    switch (opcode)
    {
    case 0:
        regs[r1] += regs[r2];
        break;
    case 1:
        regs[r1] -= regs[r2];
        break;
    case 2:
        regs[r1] *= regs[r2];
        break;
    case 3:
        if (regs[r2])
            regs[r1] /= regs[r2];
        break;
    case 4:
        if (regs[r2])
            regs[r1] %= regs[r2];
        break;
    case 5:
        regs[r1]++;
        break;
    case 6:
        regs[r1]--;
        break;
    case 7:
        regs[r1] = regs[r2];
        break;
    case 8:
        regs[r1] = ~regs[r1];
        break;
    case 9:
        regs[r1] = (regs[r1] > regs[r2]) ? 1 : 0;
        break;
    case 10:
        regs[r1] = (regs[r1] < regs[r2]) ? 1 : 0;
        break;
    case 11:
        regs[r1] = (regs[r1] == regs[r2]) ? 1 : 0;
        break;
    case 12:
        return regs[r1] % NUM_INSTRUCOES; // JMP
    case 13:
        return (regs[r1] == 0) ? regs[r2] % NUM_INSTRUCOES : pc + 1; // JZ
    case 14:
        return (regs[r1] != 0) ? regs[r2] % NUM_INSTRUCOES : pc + 1; // JNZ
    case 15:
        break; // NOP
    }
    return pc + 1;
}

int avaliar(Individuo *ind)
{
    int regs[NUM_REGS];
    for (int i = 0; i < NUM_REGS; i++)
    {
        int pos = NUM_INSTRUCOES * BITS_INSTRUCAO + i * BITS_REG;
        int local = 0;
        for (int b = 0; b < BITS_REG; b++)
        {
            local = (local << 1) | ind->genes[pos + b];
        }
        regs[i] = entrada_global[i] + local;
    }
    int pc = 0, steps = 0;
    while (steps++ < MAX_STEPS && pc < NUM_INSTRUCOES)
    {
        pc = executa_instrucao(&ind->genes[pc * BITS_INSTRUCAO], regs, pc);
    }
    if (steps >= MAX_STEPS)
        return PENALIDADE_LOOP;
    return calcula_fitness_F1(regs);
}

void salvar_todas_geracoes(int geracao)
{
    char nome[64];
    snprintf(nome, sizeof(nome), "geracao_%d.txt", geracao);
    FILE *f = fopen(nome, "w");
    const char *ops[] = {
        "ADD", "SUB", "MUL", "DIV", "MOD", "INC", "DEC", "MOV",
        "NOT", "GT", "LT", "EQ", "JMP", "JZ", "JNZ", "NOP"};

    for (int idx = 0; idx < TAMPOP; idx++)
    {
        Individuo *ind = &pop[idx];
        int regs[NUM_REGS];
        int locais[NUM_REGS];

        fprintf(f, "Individuo %d | Fitness: %d\n", idx, ind->fitness);

        fprintf(f, "Registradores globais: ");
        for (int i = 0; i < NUM_REGS; i++)
            fprintf(f, "R%d=%d ", i, entrada_global[i]);
        fprintf(f, "\n");

        fprintf(f, "Registradores locais (genes): ");
        for (int i = 0; i < NUM_REGS; i++)
        {
            int pos = NUM_INSTRUCOES * BITS_INSTRUCAO + i * BITS_REG;
            int local = 0;
            for (int b = 0; b < BITS_REG; b++)
            {
                local = (local << 1) | ind->genes[pos + b];
            }
            locais[i] = local;
            fprintf(f, "R%d=%d ", i, local);
        }
        fprintf(f, "\n");

        fprintf(f, "Registradores iniciais: ");
        for (int i = 0; i < NUM_REGS; i++)
        {
            regs[i] = entrada_global[i] + locais[i];
            fprintf(f, "R%d=%d ", i, regs[i]);
        }
        fprintf(f, "\n");

        fprintf(f, "Instrucoes:\n");
        int pc = 0, steps = 0;
        while (steps++ < NUM_INSTRUCOES * 2 && pc < NUM_INSTRUCOES)
        {
            int pos = pc * BITS_INSTRUCAO;
            int opcode = (ind->genes[pos] << 3) | (ind->genes[pos + 1] << 2) |
                         (ind->genes[pos + 2] << 1) | ind->genes[pos + 3];
            int r1 = (ind->genes[pos + 4] << 1) | ind->genes[pos + 5];
            int r2 = (ind->genes[pos + 6] << 1) | ind->genes[pos + 7];
            opcode %= 16;
            fprintf(f, "%s R%d, R%d\n", ops[opcode], r1, r2);
            pc = executa_instrucao(&ind->genes[pos], regs, pc);
        }

        fprintf(f, "Registradores finais: ");
        for (int i = 0; i < NUM_REGS; i++)
            fprintf(f, "R%d=%d ", i, regs[i]);
        fprintf(f, "\n--------------------------\n");
    }

    fclose(f);
}

int torneio()
{
    int a = rand() % TAMPOP;
    int b = rand() % TAMPOP;
    return (pop[a].fitness < pop[b].fitness) ? a : b;
}

void avaliar_populacao()
{
    for (int i = 0; i < TAMPOP; i++)
    {
        pop[i].fitness = avaliar(&pop[i]);
    }
}

void copiar_individuo(Individuo *src, Individuo *dst)
{
    memcpy(dst->genes, src->genes, CHROMOSOME_LENGTH);
    dst->fitness = src->fitness;
}

void nova_geracao()
{
    int i = 0;
    int melhor = 0;
    for (int j = 1; j < TAMPOP; j++)
    {
        if (pop[j].fitness <= pop[melhor].fitness)
            melhor = j;
    }
    int regs_temp[NUM_REGS];
    for (int i = 0; i < NUM_REGS; i++)
        regs_temp[i] = entrada_global[i];
    int pc = 0, steps = 0;
    while (steps++ < NUM_INSTRUCOES * 2 && pc < NUM_INSTRUCOES)
        pc = executa_instrucao(&pop[melhor].genes[pc * BITS_INSTRUCAO], regs_temp, pc);
    for (int i = 0; i < NUM_REGS; i++)
        entrada_global[i] = regs_temp[i];

    if (ELITISMO)
        copiar_individuo(&pop[melhor], &nova_pop[i++]);
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
        nova_pop[i++] = f1;
        if (i < TAMPOP)
            nova_pop[i++] = f2;
    }
    for (int i = 0; i < TAMPOP; i++)
        copiar_individuo(&nova_pop[i], &pop[i]);
}

int main()
{
    srand(time(NULL));
    for (int i = 0; i < NUM_REGS; i++)
        entrada_global[i] = rand() % (int)pow(2, BITS_REG);
    ;

    for (int i = 0; i < TAMPOP; i++)
        gera_individuo(&pop[i]);
    avaliar_populacao();

    for (int g = 0; g < MAXGER; g++)
    {
        salvar_todas_geracoes(g);
        nova_geracao();
        avaliar_populacao();
        int melhor = 0;
        for (int i = 1; i < TAMPOP; i++)
        {
            if (pop[i].fitness < pop[melhor].fitness)
                melhor = i;
        }
        printf("Ger %d | Melhor fitness: %d\n", g, pop[melhor].fitness);
    }

    return 0;
}
