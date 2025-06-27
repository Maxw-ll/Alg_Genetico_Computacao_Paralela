//Universidade Federal do Maranhão
//Computação Paralela
//Alunos: André Luiz Ribeiro de Araujo Lima e Mawxell Pires Silva


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define NUM_REGS 4          //numero de registradores locais
#define NUM_INSTRUCOES 4   //numero de instruções por individuo
#define BITS_INSTRUCAO 8    //quantidade de bits por instrução, sendo 4 para operação, e 4 para os dois operandos
#define BITS_REG 16          //quantidade de bits de cada registrador
#define CHROMOSOME_LENGTH (NUM_INSTRUCOES * BITS_INSTRUCAO + NUM_REGS * BITS_REG)   //tamanho de cada individuo
#define TAMPOP 30           //tamanho da populaçao
#define MAXGER 5            //numero de gerações
#define TAXCRUZ 0.8         //taxa de crossover
#define TAXMUT 0.05         //taxa de mutação
#define ELITISMO 1          //quantidade de elitismo

/*
Tabela de Operações:
Codigo | Instrução | Descrição
--------------------------------------
  0  | ADD R1,R2   | R1 = R1 + R2
  1  | SUB R1,R2   | R1 = R1 - R2
  2  | MUL R1,R2   | R1 = R1 * R2
  3  | DIV R1,R2   | R1 = R1 / R2 (if R2 != 0)
  4  | MOD R1,R2   | R1 = R1 % R2 (if R2 != 0)
  5  | INC R1      | R1 = R1 + 1
  6  | DEC R1      | R1 = R1 - 1
  7  | MOV R1,R2   | R1 = R2
  8  | AND R1,R2   | R1 = R1 & R2
  9  | OR R1,R2    | R1 = R1 | R2
 10  | XOR R1,R2   | R1 = R1 ^ R2
 11  | NOT R1      | R1 = ~R1
 12  | GT R1,R2    | R1 = (R1 > R2)
 13  | LT R1,R2    | R1 = (R1 < R2)
 14  | EQ R1,R2    | R1 = (R1 == R2)
 15  | NADA        | ------------
*/


typedef struct {
    unsigned char genes[CHROMOSOME_LENGTH];
    int aptidao;
} Individuo;

Individuo pop[TAMPOP], nova_pop[TAMPOP];
int entrada_global[NUM_REGS];

void gera_individuo(Individuo *ind) {
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        ind->genes[i] = rand() % 2;
    }
    ind->aptidao = 0;
}

void mutacao(Individuo *ind) {
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        if (((float)rand() / RAND_MAX) < TAXMUT)
            ind->genes[i] ^= 1;
    }
}

void cruzamento(Individuo pai1, Individuo pai2, Individuo *filho1, Individuo *filho2) {
    int ponto = rand() % CHROMOSOME_LENGTH;
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        if (i < ponto) {
            filho1->genes[i] = pai1.genes[i];
            filho2->genes[i] = pai2.genes[i];
        } else {
            filho1->genes[i] = pai2.genes[i];
            filho2->genes[i] = pai1.genes[i];
        }
    }
}

int calcula_fitness_F1(int *regs) {
    int resultado = regs[0] + regs[1] - regs[2] - regs[3];
    return abs(resultado);
}

void executa_instrucao(unsigned char *instrucao, int *regs) {
    int opcode = (instrucao[0]<<3) | (instrucao[1]<<2) | (instrucao[2]<<1) | instrucao[3];
    int r1 = (instrucao[4]<<1) | instrucao[5];
    int r2 = (instrucao[6]<<1) | instrucao[7];

    switch (opcode % 16) {
        case 0: regs[r1] += regs[r2]; break;
        case 1: regs[r1] -= regs[r2]; break;
        case 2: regs[r1] *= regs[r2]; break;
        case 3: if (regs[r2] != 0) regs[r1] /= regs[r2]; else regs[r1] += 7; break;
        case 4: if (regs[r2] != 0) regs[r1] %= regs[r2]; else regs[r1] += 5; break;
        case 5: regs[r1]++; break;
        case 6: regs[r1]--; break;
        case 7: regs[r1] = regs[r2]; break;
        case 8: regs[r1] = regs[r1] & regs[r2]; break;
        case 9: regs[r1] = regs[r1] | regs[r2]; break;
        case 10: regs[r1] = regs[r1] ^ regs[r2]; break;
        case 11: regs[r1] = ~regs[r1]; break;
        case 12: regs[r1] = (regs[r1] > regs[r2]) ? 1 : 0; break;
        case 13: regs[r1] = (regs[r1] < regs[r2]) ? 1 : 0; break;
        case 14: regs[r1] = (regs[r1] == regs[r2]) ? 1 : 0; break;
        case 15: // Não faça Nada 
        break;
    }
    for (int j = 0; j < NUM_REGS; j++) regs[j] &= 0x07;
}

int avaliar(Individuo *ind) {
    int regs[NUM_REGS];
    for (int i = 0; i < NUM_REGS; i++) {
        int pos = NUM_INSTRUCOES * BITS_INSTRUCAO + i * BITS_REG;
        int local = (ind->genes[pos]<<2) | (ind->genes[pos+1]<<1) | ind->genes[pos+2];
        regs[i] = (entrada_global[i] + local) & 0x07;
    }
    for (int i = 0; i < NUM_INSTRUCOES; i++) {
        executa_instrucao(&ind->genes[i * BITS_INSTRUCAO], regs);
    }
    return calcula_fitness_F1(regs);
}



int torneio() {
    int a = rand() % TAMPOP;
    int b = rand() % TAMPOP;
    return (pop[a].aptidao < pop[b].aptidao) ? a : b;
}

void avaliar_populacao() {
    for (int i = 0; i < TAMPOP; i++) {
        pop[i].aptidao = avaliar(&pop[i]);
    }
}

void copiar_individuo(Individuo *src, Individuo *dst) {
    memcpy(dst->genes, src->genes, CHROMOSOME_LENGTH);
    dst->aptidao = src->aptidao;
}

void salvar_todas_geracoes(int geracao) {
    char nome[64];
    snprintf(nome, sizeof(nome), "geracao_%d.txt", geracao);
    FILE *f = fopen(nome, "w");
    const char *ops[] = {
        "ADD", "SUB", "MUL", "DIV", "MOD", "INC", "DEC", "MOV",
        "AND", "OR", "XOR", "NOT", "GT", "LT", "EQ", "NE"
    };

    for (int idx = 0; idx < TAMPOP; idx++) {
        Individuo *ind = &pop[idx];
        int regs[NUM_REGS];
        int locais[NUM_REGS];

        fprintf(f, "Individuo %d | Fitness: %d\n", idx, ind->aptidao);

        fprintf(f, "Registradores globais: ");
        for (int i = 0; i < NUM_REGS; i++)
            fprintf(f, "R%d=%d ", i, entrada_global[i]);
        fprintf(f, "\n");

        fprintf(f, "Registradores locais (genes): ");
        for (int i = 0; i < NUM_REGS; i++) {
            int pos = NUM_INSTRUCOES * BITS_INSTRUCAO + i * BITS_REG;
            locais[i] = (ind->genes[pos] << 2) | (ind->genes[pos + 1] << 1) | ind->genes[pos + 2];
            fprintf(f, "R%d=%d ", i, locais[i]);
        }
        fprintf(f, "\n");

        fprintf(f, "Registradores iniciais: ");
        for (int i = 0; i < NUM_REGS; i++) {
            regs[i] = (entrada_global[i] + locais[i]) & 0x07;
            fprintf(f, "R%d=%d ", i, regs[i]);
        }
        fprintf(f, "\n");

        fprintf(f, "Instrucoes:\n");
        for (int i = 0; i < NUM_INSTRUCOES; i++) {
            int pos = i * BITS_INSTRUCAO;
            int opcode = (ind->genes[pos] << 3) | (ind->genes[pos + 1] << 2) |
                         (ind->genes[pos + 2] << 1) | ind->genes[pos + 3];
            int r1 = (ind->genes[pos + 4] << 1) | ind->genes[pos + 5];
            int r2 = (ind->genes[pos + 6] << 1) | ind->genes[pos + 7];
            fprintf(f, "%s R%d, R%d\n", ops[opcode % 16], r1, r2);
            executa_instrucao(&ind->genes[pos], regs);
        }

        fprintf(f, "Registradores finais: ");
        for (int i = 0; i < NUM_REGS; i++)
            fprintf(f, "R%d=%d ", i, regs[i]);
        fprintf(f, "\n--------------------------\n");
    }

    fclose(f);
}



void nova_geracao() {
    int i = 0;
    int melhor = 0;
    for (int j = 1; j < TAMPOP; j++) {
        if (pop[j].aptidao <= pop[melhor].aptidao) melhor = j;
    }
    int regs_temp[NUM_REGS];
    for (int i = 0; i < NUM_REGS; i++) regs_temp[i] = entrada_global[i];
    for (int i = 0; i < NUM_INSTRUCOES; i++) {
        executa_instrucao(&pop[melhor].genes[i * BITS_INSTRUCAO], regs_temp);
    }
    for (int i = 0; i < NUM_REGS; i++) entrada_global[i] = regs_temp[i];

    if (ELITISMO) copiar_individuo(&pop[melhor], &nova_pop[i++]);
    while (i < TAMPOP) {
        int p1 = torneio();
        int p2 = torneio();
        Individuo f1, f2;
        if (((float)rand() / RAND_MAX) < TAXCRUZ) {
            cruzamento(pop[p1], pop[p2], &f1, &f2);
        } else {
            copiar_individuo(&pop[p1], &f1);
            copiar_individuo(&pop[p2], &f2);
        }
        mutacao(&f1);
        mutacao(&f2);
        nova_pop[i++] = f1;
        if (i < TAMPOP) nova_pop[i++] = f2;
    }
    for (int i = 0; i < TAMPOP; i++) copiar_individuo(&nova_pop[i], &pop[i]);
}

int main() {
    srand(time(NULL));
    for (int i = 0; i < NUM_REGS; i++) entrada_global[i] = rand() % 65536;

    for (int i = 0; i < TAMPOP; i++) gera_individuo(&pop[i]);
    avaliar_populacao();

    for (int g = 0; g < MAXGER; g++) {
        salvar_todas_geracoes(g);
        nova_geracao();
        avaliar_populacao();
        int melhor = 0;
        for (int i = 1; i < TAMPOP; i++) {
            if (pop[i].aptidao < pop[melhor].aptidao) melhor = i;
        }
        printf("Ger %d | Melhor fitness: %d\n", g, pop[melhor].aptidao);
    }

    return 0;
}
