#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

/* Data strucutres */
// polinomios 
typedef struct {
        int power;
        double * coef;
} poly_t, *poly;

// chain de polinomios
struct chain_t {
    poly pol;
    struct chain_t *seg;
};
typedef struct chain_t chain;

// pares
struct delta_t{
    double point;
    int count;
};
typedef struct delta_t delta;

// chain de pares
struct interval_t{
    delta pair;
    struct interval_t *seg;
};
typedef struct interval_t interval;

// tuplas
struct tupla_t{
    double a;
    double b;
    struct tupla_t *seg;
};
typedef struct tupla_t tupla;

// lista de double
struct val_t {
    double val;
    struct val_t *seg;
};
typedef struct val_t val;

/* Definicao */
#define E(x, i) (x)->coef[i]

/* Metodos especificos para C  */ 
 
/* passing in negative power to have a zeroed poly */
poly p_new(int power, ...)
{
        int i, zeroed = 0;
        va_list ap;
 
        if (power < 0) {
                power = -power;
                zeroed = 1;
        }
 
        poly p = malloc(sizeof(poly_t));
        p->power = power;
        p->coef = malloc(sizeof(double) * ++power);
 
        if (zeroed)
                for (i = 0; i < power; i++) p->coef[i] = 0;
        else {
                va_start(ap, power);
                for (i = 0; i < power; i++)
                        E(p, i) = va_arg(ap, double);
                va_end(ap);
        }
 
        return p;
}
 
void p_del(poly p)
{
        free(p->coef);
        free(p);
}
 
void p_print(poly p)
{
        int i;
        for (i = 0; i <= p->power; i++)
                printf("%g ", E(p, i));
        printf("\n");
}
 
poly p_copy(poly p)
{
        poly q = p_new(-p->power);
        memcpy(q->coef, p->coef, sizeof(double) * (1 + p->power));
        return q;
}

void chain_print(chain *sturm)
{
    chain *p;
    for(p=sturm; p != NULL;p = p->seg)
        p_print(p->pol);
}

void chain_insert(poly p, chain *c)
{
    chain *nova;
    nova = malloc(sizeof(chain));
    nova->pol=p;
    nova->seg = c->seg;
    c->seg=nova;
}

int chain_size(chain *sturm)
{
    int i = 0;
    chain *p;
    for(p=sturm; p != NULL;p = p->seg)
        i += 1;
    return i;
}

void interval_insert(delta d, interval *it){
    interval *nova;
    nova = malloc(sizeof(interval));
    nova->pair=d;
    nova->seg=it->seg;
    it->seg=nova;
}

void interval_print(interval *it)
{
    interval *p;
    for(p=it; p != NULL;p = p->seg)
        printf("interval [%.2f %d]\n", p->pair.point, p->pair.count);
}

void tupla_insert(double a, double b, tupla *it){
    tupla *nova;
    nova = malloc(sizeof(tupla));
    nova->a=a;
    nova->b=b;
    nova->seg=it->seg;
    it->seg=nova;
}

void tupla_print(tupla *it)
{
    tupla *p;
    for(p=it; p != NULL;p = p->seg)
        printf("tupla [%.2f %.2f]\n", p->a, p->b);
}


double max(int len, double *tmp)
{
    double maxValue=fabs(tmp[0]/tmp[len-1]);
    for(int i=0; i < len-1; i++){
        if(fabs(tmp[i]/tmp[len-1]) < fabs(tmp[i+1]/tmp[len-1])){
            maxValue = fabs(tmp[i+1]/tmp[len-1]);
        }
    }
    return maxValue;
}

void val_print(val it)
{
    val *p;
    for(p=&it; p != NULL;p = p->seg)
        printf("%.5f ", p->val);
}

void val_insert(double a, val *it){
    val *nova;
    nova = malloc(sizeof(val));
    nova->val=a;
    nova->seg=it->seg;
    it->seg=nova;
}

/* Metodos convertidos  */ 

/* p: poly;  d: divisor;  r: remainder; returns quotient */
poly p_div(poly p, poly d, poly* r)
{

    poly q;
    int i, j;
    int power = p->power - d->power;
    double ratio;

    if (power < 0){
        q = p_new(-1);
        *r = p_copy(q);
        return q;
    }

    q = p_new(-power);
    *r= p_copy(p);

    for (i = p->power; i >= d->power; i--) {
            E(q, i - d->power) = ratio = E(*r, i) / E(d, d->power);
            E(*r ,i) = 0;

            for (j = 0; j < d->power; j++)
                    E(*r, i - d->power + j) -= E(d, j) * ratio;
    }
    while (! E(*r, --(*r)->power));

    return q;

}

poly p_dev(poly p)
{
    int i;
    int power = p->power - 1;
    if(power < 0) return 0;

    poly q = p_new(-power);
    for(i = q->power; i >= 0; i--){
        E(q, i)=E(p, i+1)*(i+1);
    }
    return q;
}

double p_eval(poly p, double x)
{
    double tmp = 0;
    int i;
    for(i=p->power; i >= 0; i--){
        if(i == 0){
            tmp += E(p, i);
        }
        else{
            tmp += E(p, i)*pow(x, (double)i);
        }
    }
    return tmp;
}

int p_countSignChanges(poly p)
{
    int i, count=0;
    for(i=p->power; i >= 1; i--){
        if((E(p, i)>=0 && E(p, i-1)<0) || (E(p, i)<0 && E(p, i-1)>=0)){    
            count += 1;
        }
    }
    return count;
}

int countSignChanges(int len, double *tmp)
{
    int i, count=0; 
    for(i=0; i < len-1; i++){
        if((tmp[i]>=0 && tmp[i+1]<0) || (tmp[i]<0 && tmp[i+1]>=0)){
                count += 1;
        }
    }
    return count;
}

chain sturm_chain(poly p)
{
    chain c;
    c.pol = p;
    c.seg = NULL;
    
    poly dv = p_dev(p);
    chain_insert(dv, &c);    
    
    chain *a, *b;
    a=&c;
    b=a->seg;
    while(b != NULL){
        if(b->pol->power > 0){
            poly remainder;
            p_div(a->pol, b->pol, &remainder);
            for (int i = 0; i <= remainder->power; i++){
              E(remainder, i) = -1*E(remainder, i);
            }
            chain_insert(remainder, b);
        }
        a=b;
        b=b->seg;
    }

    return c;
}


tupla verify_root(interval *it)
{
    int change;
    interval *p;

    tupla ac;
    ac.a=0.0;
    ac.b=0.0;
    ac.seg=NULL;

    tupla *acc;
    acc=&ac;

    for(p=it;p != NULL;p = p->seg){
        if(p->seg != NULL){
            change = p->pair.count - p->seg->pair.count;
            if( change > 0 ){
                tupla_insert(p->pair.point, p->seg->pair.point, acc);
                acc=acc->seg;
            }
        }
    }
    return *ac.seg;
}


tupla isolate_all_roots(poly polinomio, double min, double max)
{
    double h = 1;
    double point = min;
    chain s_chain = sturm_chain(polinomio);
        
    int len = chain_size(&s_chain);
    
    chain *p;
    int i=0;

    double *tmp = malloc(len * sizeof(double));

    for(p=&s_chain; p != NULL;p = p->seg){
        tmp[i] = p_eval(p->pol, point);
        i++;
    }

    delta d;
    d.point=point;
    d.count=countSignChanges(len, tmp);
    free(tmp);

    interval it;
    it.pair = d;
    it.seg = NULL;
    
    point += h;

    interval *it_a;
    it_a=&it;

    while(point <= max){
        i=0;
        double *tmp = malloc(len * sizeof(double));

        for(p=&s_chain; p != NULL;p = p->seg){
            tmp[i] = p_eval(p->pol, point);
            i++;
        }

        delta d;
        d.point=point;
        d.count=countSignChanges(i, tmp);
        free(tmp); 

        interval_insert(d, it_a);
        it_a=it_a->seg;
        
        point += h;
    }
    
    tupla bar = verify_root(&it);
    
    return bar;
}

double newton_raphson(poly p, tupla t)
{
    double x = (t.b - t.a)/2 + t.a;

    int n = 100;
    double next_ite;
    double accuracy=0.0001;

    for(int i=0;i < n; i++){
        if(p_eval(p_dev(p), x) == 0)
            break;
        else{
            next_ite = x - (p_eval(p, x)/p_eval(p_dev(p), x));
            if(fabs(x-next_ite) < accuracy && fabs(p_eval(p, next_ite)) < accuracy){
                x=next_ite;
                break;
            }
            x=next_ite;
        }
    }

    //printf("[!] %.5f\n", p_eval(p, x));
    return x;
}


double root_radius(poly p)
{
    double resp = 3.0 + max(p->power+1, p->coef);
    return resp;
}


val find_all_roots(poly p)
{
    chain s_chain2 = sturm_chain(p);
    //printf("        poly: "); p_print(p);
    //chain_print(&s_chain2);

    double r = root_radius(p);

    tupla foo3 = isolate_all_roots(p, -r, r);
    //tupla_print(&foo3);

    tupla *t;
    
    val root;
    root.val = 0.0;
    root.seg = NULL;
    val *roots;
    roots=&root;
    
    for(t=&foo3; t != NULL;t = t->seg){
        val_insert(newton_raphson(p, *t), roots);
    }

    p_del(p); 
    
    return *(roots->seg);
}

int main(int argc, char **argv)
{
    FILE *f = fopen(argv[1], "r"); // "r" for read
    
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    
    val root;
    root.val = 0.0;
    root.seg = NULL;
    val *roots;
    roots=&root;
    
    int lines=0;
    while ((read = getline(&line, &len, f)) != -1) {
        printf("%s", line);
        double valor;
        sscanf(line, "%lf", &valor);
        printf("valor %.2f\n", valor);
        lines++;
        val_insert(valor, roots);
    }
    fclose(f);  /* close the file */ 
    printf("[lines] %d\n", lines);
    
    val *p;
    int i = lines-1;
    poly polinomio = p_new(-i);
    p_print(polinomio);
    for(p=roots->seg; p != NULL;p = p->seg){
        //printf("%.5f ", p->val);
        E(polinomio, i) = p->val;
        i--;
    }
    p_print(polinomio);
    printf("%d", polinomio->power);
    
/*
    printf("[*] TESTES\n");

    poly p1 = p_new(3, -24., 6., -4., 1.);
    poly d1 = p_new(2, -1., 0., 1.);
    poly r1;
    poly q1 = p_div(p1, d1, &r1);
    printf("\n[poldiv]\n");
    printf("[-24,6,-4,1]    poly: "); p_print(p1);
    printf("[-1,0,1]        div:  "); p_print(d1);
    printf("[-4,1]          quot: "); p_print(q1);
    printf("[-28, 7]        rem:  "); p_print(r1);    
    p_del(p1);
    p_del(q1);
    p_del(r1);
    p_del(d1);

    poly p2 = p_new(5, -1., -3., 0., 0., 0., 1.);
    printf("\n[sturm_chain]\n");
    chain s_chain1 = sturm_chain(p2);
    printf("        poly: "); p_print(p2);
    chain_print(&s_chain1);
    p_del(p2);

    printf("\n[countSignChanges]\n");
    poly p3 = p_new(1, 0., -1.);
    poly q3 = p_new(2, -1., 0., -1.);
    printf("        poly: "); p_print(p3);
    printf("        count: %d\n", p_countSignChanges(p3));
    printf("        poly: "); p_print(q3);
    printf("        count: %d\n", p_countSignChanges(q3));
    p_del(p3);
    p_del(q3);    

    printf("\n[isolate_all_roots]\n");
    poly p4 = p_new(2, -4., 0., 1.);
    tupla foo1 = isolate_all_roots(p4, -10, 10);
    tupla_print(&foo1);
    p_del(p4);

    printf("\n[isolate_all_roots #2]\n");
    poly p5 = p_new(16, 1.,2.,3.,4.,56.,7.,8.,8.,9.,6.,5.,4.,2.,2.,42.,423.,3.);
    tupla foo2 = isolate_all_roots(p5, -150, 150);
    tupla_print(&foo2);
    p_del(p5);
*/
    printf("\n[find_all_roots]\n");
    poly final = p_new(16, 1.,2.,3.,4.,56.,7.,8.,8.,9.,6.,5.,4.,2.,2.,42.,423.,3.);
    printf("[polinomio]: ");p_print(final);
    val t = find_all_roots(final);
    printf("[raizes]: ");val_print(t);
    printf("\n");

    return 0;
}
