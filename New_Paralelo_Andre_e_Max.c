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
#define TAMPOP 10
/*Quantidade Máxima de Gerações*/
#define MAXGER 20
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

/*Estrutura do Indivíduo Composto pelo Cromossomo (A Parte que é evoluída) e seus Registradores, além de sua aptidão*/
typedef struct
{
    unsigned char genes[CHROMOSOME_LENGTH];
    int regs[NUM_REGS];
    int age_regs[NUM_REGS]; // Apenas para fim de anotação dos valores antigos e atualizados de registradores
    int age_pcs[MAX_STEPS]; // Armazenar o PC apos executar, apenas para fim de anotação no txt
    int fitness;
} Individuo;

/*Inicia uma população de Indivíduos e Nova, que servirá para adicionar os indívudos novos com base na anterior. Depois nova_pop é copiada para pop*/
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
        {
            /*Sorteia uma Instrução e reescreve de forma aleatória*/
            int instr_index = rand() % NUM_INSTRUCOES;
            int base = instr_index * BITS_INSTRUCAO;
            for (int j = 0; j < BITS_INSTRUCAO; j++)
            {
                ind->genes[base + j] = rand() % 2;
            }
        }
        break;
        case 2:
        {
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
    case F1: /* A + B = C + D*/
        resultado = abs((regs[0] + regs[1]) - (regs[2] + regs[3]));
        break;
    case F2: /* A + B > C - D*/
        if ((regs[0] + regs[1]) > regs[2] - regs[3])
        {
            resultado = 0;
        }
        else
        {
            resultado = (int)pow(2, BITS_REG) - 1;
        }
        break;
    case F3: /* A % B == 0 -> B divide A*/
        if (regs[1] != 0 && regs[0] % regs[1] == 0)
        {
            resultado = 0;
        }
        else
        {
            resultado = (int)pow(2, BITS_REG) - 1;
        }
        break;
    case F4: /* A=B+1 && B=C+1 && C=D+1 */
        if (regs[0] == regs[1] + 1 && regs[1] == regs[2] + 1 && regs[2] == regs[3] + 1)
        {
            resultado = 0;
        }
        else
        {
            resultado = (int)pow(2, BITS_REG) - 1;
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
        break;
    }
    }
    return resultado;
}

/*Função para Normalizar os valores dos registradores, uma vez que como cada registrador é um inteiro em si, precisamos garantir que o valor dentro dele esteja no intervalo correto*/
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

/*
Tabela de Operações:
Código | Instrução   | Descrição
------------------------------------------
  0    | ADD R1,R2   | R1 = R1 + R2
  1    | SUB R1,R2   | R1 = R1 - R2
  2    | MUL R1,R2   | R1 = R1 * R2
  3    | DIV R1,R2   | R1 = R1 / R2 (se R2 ≠ 0)
  4    | MOD R1,R2   | R1 = R1 % R2 (se R2 ≠ 0)
  5    | INC R1      | R1 = R1 + 1
  6    | DEC R1      | R1 = R1 - 1
  7    | MOV R1,R2   | R1 = R2
  8    | NOT R1      | R1 = ~R1
  9    | GT R1,R2    | R1 = (R1 > R2) ? 1 : 0
 10    | LT R1,R2    | R1 = (R1 < R2) ? 1 : 0
 11    | EQ R1,R2    | R1 = (R1 == R2) ? 1 : 0
 12    | JMP R1      | Salta para instrução no índice R1 % NUM_INSTRUCOES
 13    | JZ R1,R2    | Se R1 == 0, salta para índice R2 % NUM_INSTRUCOES; senão, PC++
 14    | JNZ R1,R2   | Se R1 ≠ 0, salta para índice R2 % NUM_INSTRUCOES; senão, PC++
 15    | NOP         | Não faz nada
*/

/*Executa uma Instrução com base na Tabela*/
int executa_instrucao(unsigned char *instrucao, int *regs, int pc)
{
    // Transformamos os bits em inteiros. Os 4 primeiros bits geram o opcode e os 4 seguintes, os registradores que serão utilizados na instrução
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
        if (regs[r2])
        {
            regs[r1] /= regs[r2];
        }
        break;
    case 4: // MOD
        if (regs[r2])
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
    case 15:
        break; // NOP
    }
    normalized_regs(regs);
    return pc + 1;
}

/*Zerar o vetor que armazena o PC*/
void zerar_pcs(Individuo *ind)
{
    for (int i = 0; i < MAX_STEPS; i++)
    {
        ind->age_pcs[i] = PENALIDADE_LOOP;
    }
}

/*Função para avaliar um Indivíduo -> Executar suas instruções e calcular o Fitness*/
int avaliar(Individuo *ind, int Tipo_Fitness)
{

    int pc = 0, steps = 0;

    zerar_pcs(ind);

    /*Copiar os valores antigos dos registradores para fim de anotação*/
    for (int i = 0; i < NUM_REGS; i++)
    {
        ind->age_regs[i] = ind->regs[i];
    }

    while (steps < MAX_STEPS && pc < NUM_INSTRUCOES)
    {
        ind->age_pcs[steps] = pc;
        pc = executa_instrucao(&ind->genes[pc * BITS_INSTRUCAO], ind->regs, pc);
        steps++;
    }

    if (steps >= MAX_STEPS)
    {
        return PENALIDADE_LOOP;
    }

    return calcula_fitness(ind->regs, Tipo_Fitness);
}


/* Ao mesmo tempo que executa as instruções de cada indivíduo, salva suas informações*/
void salvar_todas_geracoes(int geracao)
{
    const char *ops[] = {
        "ADD", "SUB", "MUL", "DIV", "MOD", "INC", "DEC", "MOV",
        "NOT", "GT", "LT", "EQ", "JMP", "JZ", "JNZ", "NOP"};

    char nome[64];
    snprintf(nome, sizeof(nome), "geracao_%d.txt", geracao);
    FILE *f = fopen(nome, "w");

    for (int idx = 0; idx < TAMPOP; idx++)
    {
        Individuo *ind = &pop[idx];

        fprintf(f, "Individuo %d \n", idx);

        fprintf(f, "Cromossomo: ");

        for (int i = 0; i < CHROMOSOME_LENGTH; i++)
        {
            fprintf(f, "%d", ind->genes[i]);
        }
        fprintf(f, "\nRegistradores -> Valores Iniciais: ");
        for (int i = 0; i < NUM_REGS; i++)
        {
            fprintf(f, "R%d=%d ", i, ind->age_regs[i]);
        }
        fprintf(f, "\n");

        fprintf(f, "Instrucoes:\n");

        int pc = 0, steps = 0;

        // if (idx == 1)
        // {
        //     for(int i=0; i<MAX_STEPS; i++)
        //     {
        //         printf("PC[%d] -> %d ", i, ind->age_pcs[i]);
        //     }
        // }
        /*Executa as instruções do indivíduo e anota o processo*/
        while (steps < MAX_STEPS && pc < NUM_INSTRUCOES)
        {
            int pos = pc * BITS_INSTRUCAO;
            int opcode = (ind->genes[pos] << 3) | (ind->genes[pos + 1] << 2) | (ind->genes[pos + 2] << 1) | ind->genes[pos + 3];
            int r1 = (ind->genes[pos + 4] << 1) | ind->genes[pos + 5];
            int r2 = (ind->genes[pos + 6] << 1) | ind->genes[pos + 7];
            
            fprintf(f, "%s R%d, R%d\n", ops[opcode], r1, r2);
            steps++;
            pc = ind->age_pcs[steps];
        }

        zerar_pcs(ind);

        fprintf(f, "Registradores -> Valores Finais: ");
        for (int i = 0; i < NUM_REGS; i++)
        {
            fprintf(f, "R%d=%d ", i, ind->regs[i]);
        }
        fprintf(f, "\n");
        fprintf(f, "Fitness Final: %d\n", ind->fitness);
        fprintf(f, "\n--------------------------\n");
    }

    fclose(f);
}

/* Torneio para Escolher melhor Indivíduo dentre dois escolhidos aleatoriamente*/
int torneio()
{
    int a = rand() % TAMPOP;
    int b = rand() % TAMPOP;
    return (pop[a].fitness < pop[b].fitness) ? a : b;
}

/*Aplicar o fitness em toda a população*/
void avaliar_populacao(int Tipo_Fitness)
{
    for (int i = 0; i < TAMPOP; i++)
    {
        pop[i].fitness = avaliar(&pop[i], Tipo_Fitness);
    }
}

/*Copiar Indivíduo para outro Novo*/
void copiar_individuo(Individuo *src, Individuo *dst)
{
    memcpy(dst->genes, src->genes, CHROMOSOME_LENGTH);
    dst->fitness = src->fitness;
}

/* Após cada Nova geração ser gerada, os valores de todos os registradores são aleatorizados*/
void gerar_novas_entradas()
{
    for (int i = 0; i < TAMPOP; i++)
    {
        for (int j = 0; j < NUM_REGS; j++)
        {
            if (rand() % 2 == 0)
            {
                pop[i].regs[j] += ((rand() % 10) + 1);
            }
            else
            {
                pop[i].regs[j] -= ((rand() % 10) + 1);
            }
        }

        normalized_regs(pop[i].regs);
    }
}

/*Gerar uma Nova geração a partir da Anterior*/
void nova_geracao()
{
    int i = 0;
    int melhor = 0;
    for (int j = 1; j < TAMPOP; j++)
    {
        if (pop[j].fitness <= pop[melhor].fitness)
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
            nova_pop[i] = f2;
        i++;
    }

    for (int i = 0; i < TAMPOP; i++)
        copiar_individuo(&nova_pop[i], &pop[i]);

    gerar_novas_entradas();
}

int main()
{
    srand(time(NULL));

    int tipo_fitness = F1;

    for (int i = 0; i < TAMPOP; i++)
    {
        gera_individuo(&pop[i]);
    }

    avaliar_populacao(tipo_fitness);

    /*Arquivo de Log*/
    char nome[64];
    snprintf(nome, sizeof(nome), "Log_Geracoes.txt");
    FILE *f = fopen(nome, "w");

    for (int g = 0; g < MAXGER; g++)
    {
        salvar_todas_geracoes(g);
        int melhor = 0;
        for (int i = 1; i < TAMPOP; i++)
        {
            if (pop[i].fitness < pop[melhor].fitness)
                melhor = i;
        }

        printf("Ger %d | Melhor fitness: %d\n", g, pop[melhor].fitness);
        fprintf(f, "Ger %d | Melhor fitness: %d\n", g, pop[melhor].fitness);

        nova_geracao();
        avaliar_populacao(tipo_fitness);
    }
    fclose(f);
    return 0;
}
