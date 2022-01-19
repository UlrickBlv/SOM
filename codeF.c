#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define ROUGE "\033[0;31m"
#define VERT "\033[0;32m"
#define VIOLET "\033[0;35m"




struct vecteur{
double *arr; /* valeur des donnes iris.data */
char *etiquette;/* etiquette avec les especes */
double norme; /* calcule des normes */
};

struct node
{
double act;/* distance euclidienne */
char *id; /* id des especes */
double *w; /* vecteur de poids*/
};
typedef struct node node;


struct bmu{
int colonne; /* colonne du bmu */
int ligne; /* ligne du bmu */
struct bmu *suivant; /* bmu suivant */
};

typedef struct bmu bmuF;


typedef struct liste_bmu{
    int nombre_bmu; /* le nombre de bmu */
    struct bmu *premier; /* premier bmu */

}liste_bmu;

struct reseau{
//double alpha;
char *etiq; /* etiquettes des especes */
double *capteur; /* vecteur courant*/
node **map; /* matrice du reseau */
int taille_voisinage;
}reseau;

struct N_config{
int taille_v; /* taille des vecteurs du fichier */
int taille_f; /* Taille du fichier */
int nombre_l; /* Nombre de ligne du reseau  */
int nombre_c; /* Nombre de colonne du reseau  */
int nb_iteration; /* Nombre d'iteration */
int nb_ite1; /* nombre d'iteration phase 1 */
double alphaM;
}N_config;

/* creation variable de la structure vecteur */
struct vecteur *tab_v;

/* declaration de variables */
double *vectM;
double *min;
double *max;
int *index_tab;



/* Normalisation des vecteurs*/
void normaliser(int i,int size)
{
    double somme=0.0;
    for(int j=0;j<size;j++)
        somme+=(tab_v[i].arr[j])*(tab_v[i].arr[j]);
    tab_v[i].norme=sqrt(somme);
}

/* Denormalisation des vecteurs*/
void denormaliser(int nb_ligne,int nb_colonne)
{
    int i,j;
    for(i=0;i<nb_ligne;i++)
    {
        for(j=0;j<nb_colonne;j++)
            tab_v[i].arr[j]/=tab_v[i].norme;
    }
}

/* Extraction des vecteurs du fichier pour les stocker */
void parser_fi(int nbLigne, int nbColonne){
    FILE *f;
    char *str=malloc(sizeof(char)*500);
    f=fopen("iris.data","r");
    for(int i=0;i<nbLigne;i++){
        fscanf(f,"%s",str);
        char *sep=strtok(str,",");
        for(int j=0;j<nbColonne;j++){
            tab_v[i].arr[j]=atof(sep);
            sep=strtok(NULL,",");
        }
    if(strcmp(sep, "Iris-setosa") == 0){
        strcpy(tab_v[i].etiquette,"Se");
    }
    else if(strcmp(sep,"Iris-versicolor") == 0)
    {
        strcpy(tab_v[i].etiquette,"Ve");
    }
    else if(strcmp(sep, "Iris-virginica") == 0)
    {
        strcpy(tab_v[i].etiquette,"Vi");
    }
    normaliser(i,nbColonne);
    
    }
    denormaliser(nbLigne,nbColonne);
    fclose(f);
    free(str);
}



/* vecteur moyen */
void vect_moyen(int nb_ligne,int nb_colonne)
{
    vectM= malloc(nb_colonne*sizeof(double));
    for(int i=0;i<nb_colonne;i++){
        for(int j=0;j<nb_ligne;j++){
            vectM[i]+=tab_v[j].arr[i];
        }
        vectM[i]/=nb_ligne;
    }
}

/* vecteur min */
void Minv(int nb_colonne){
   min=malloc(nb_colonne*sizeof(double));
   for(int i=0;i<nb_colonne;i++){
       min[i]=vectM[i]- 0.05;
    }
  }

/* vecteur max */
void Maxv(int nb_colonne){
   max=malloc(nb_colonne*sizeof(double));
   for(int i=0;i<nb_colonne;i++){
       max[i]=vectM[i]+ 0.05;
    }
}
  

double* init_alea()
{
    int i;
    double alea=(double)rand()/RAND_MAX;
    double *tmp=malloc(4*sizeof(double));

    for(i=0;i<4;i++)
        {
            tmp[i]=alea*(max[i]-min[i])+min[i];
        }
    return tmp;
}


void alloue_tab(int n)
{
    tab_v=malloc(n*sizeof(struct vecteur));
    for(int i=0;i<n;i++)
    {
        tab_v[i].arr=malloc(4*sizeof(double));
        tab_v[i].etiquette=malloc(20*sizeof(char));
    }
}

/* initialisation du tableau */
void init_tab(int n)
{
    index_tab=malloc(sizeof(int)*n);
    for(int i=0;i<n;i++)
        index_tab[i]=i;
}

/* Mélange des valeurs du tableau initialisé */
void shuffle(int Nb)
{       
    int tmp;
    int random;
    srand(time(NULL));
    for (int i=0;i<Nb;i++){
        random = rand()%Nb;   
        tmp = index_tab[i];
        index_tab[i]=index_tab[random];
        index_tab[random]=tmp;
    }
}

/* distance euclidienne */
double dts_eucl(double *vect1,double *vect2,int nb_colonne){
    double somme = 0.0;
    for(int i =0;i<nb_colonne;i++){
        somme += somme+pow((vect1[i]-vect2[i]),2);
    }
    somme = sqrt(somme);
    return somme;
}

/* Calcul du voisinage et permets de mettre à jour le Bmu  */
void voisinage(bmuF* b_mu,double alpha)
{
    int voisin=reseau.taille_voisinage;
    int X_min,X_max,Y_min,Y_max;

    for(;voisin>=0;voisin--)
    {
        if(b_mu->ligne-voisin<0)
            X_min=0;
        else
            X_min=b_mu->ligne-voisin;
        if(b_mu->colonne-voisin<0)
            Y_min=0;
        else
            Y_min=b_mu->colonne-voisin;

        if(b_mu->ligne+voisin>N_config.nombre_l-1)
            X_max=N_config.nombre_l-1;
        else
            X_max=b_mu->ligne+voisin;
        if(b_mu->colonne+voisin>N_config.nombre_c-1)
            Y_max=N_config.nombre_c-1;
        else
            Y_max=b_mu->colonne+voisin;

        for(int i=X_min;i<=X_max;i++)
            for(int j=Y_min;j<=Y_max;j++)
            {

                for(int k=0;k<4;k++)
                    {

                        reseau.map[i][j].w[k]+=alpha*(reseau.capteur[k]-reseau.map[i][j].w[k]);
                    }
            }
    }
}


/* Initialiser le réseau de neurone */
void config_NeuronMap(){
    reseau.map=malloc(N_config.nombre_l*sizeof(node*));
     for(int i=0;i<N_config.nombre_l;i++){
          reseau.map[i]=malloc(N_config.nombre_c*sizeof(node));
     }

    for(int i=0;i<N_config.nombre_l;i++){
        for(int j=0;j<N_config.nombre_c;j++){

            reseau.map[i][j].w=(double *)malloc(sizeof(double)*N_config.taille_v);
            reseau.map[i][j].w=init_alea();
            reseau.map[i][j].id=malloc(20*sizeof(char));
            strcpy(reseau.map[i][j].id, "* ");
        }
    }
}

/* Affichage de la map */
void affiche_map(){
    printf(VIOLET "Iris-Virginica\t");
    printf(VERT "Iris-Virginica\t");
    printf(ROUGE "Iris-Setosa\t \n");
    printf("\033[0m");
    printf("\n");
    for(int i=0;i<6;i++){
        for(int j=0;j<10;j++){
            if (!strcmp(reseau.map[i][j].id, "Vi"))
            {
                printf(VIOLET " " );
            }

            if (!strcmp(reseau.map[i][j].id, "Ve"))
            {
                printf(VERT " ");
            }

            if (!strcmp(reseau.map[i][j].id, "Se"))
            {
                printf(ROUGE " ");
            }
            printf("%s",reseau.map[i][j].id);
            printf("\033[0m");
        }
        printf("\n");
        
    } 
}

/* Initialisation de la liste */
liste_bmu *initialisation(){
     liste_bmu *listeBMU=malloc(sizeof(*listeBMU));
     bmuF *bmu=malloc(sizeof(*bmu));
    if(listeBMU == NULL || bmu == NULL){
        exit(EXIT_FAILURE);
    }
    bmu->colonne=0;
    bmu->ligne=0;
    bmu->suivant=NULL;
    listeBMU->premier=bmu;
    listeBMU->nombre_bmu+=1;

    return listeBMU;
}

/* Fonction insertion qui va permettre de remplir la liste */
void insertion(liste_bmu *listeBMU,int newLigne,int newColonne){
    bmuF *nouveau=malloc(sizeof(*nouveau));
    if(listeBMU == NULL || nouveau == NULL){
        exit(EXIT_FAILURE);
    }
    nouveau->colonne = newColonne;
    nouveau->ligne = newLigne;

    nouveau->suivant = listeBMU->premier;
    listeBMU->premier = nouveau;
    listeBMU->nombre_bmu+=1;

}

/* Fonction insertion qui va permettre de supprimer la liste  */
void suppression(liste_bmu *listeBMU){
    if(listeBMU == NULL){
        exit(EXIT_FAILURE);
    }
    while(listeBMU->premier->suivant){
        bmuF *aSupprimer = listeBMU->premier;
        listeBMU->premier = listeBMU->premier->suivant;
        free(aSupprimer);
    }
    listeBMU->premier=NULL;
    listeBMU->nombre_bmu=0;
}


/* Fonction qui va permettre de detecter le bmu  */
bmuF *detecte_BMU(int nb_ligne,int nb_colonne){
    int random;
    double dist;
  liste_bmu *listeBMU = initialisation();
   double dist_euc_min = 10;
    for(int i=0; i<nb_ligne;i++){
        for(int j=0; j<nb_colonne;j++){
            dist=dts_eucl(reseau.capteur,reseau.map[i][j].w,4);
            if(dist < dist_euc_min){
                //supression des valeurs de la liste
                suppression(listeBMU);
                //ajout des valeurs dans la liste
                insertion(listeBMU,i,j);
                  //initialisation nouveau distance minimun
                dist_euc_min =  dist;

            }
            else if(dist_euc_min == dist){
               //ajout des valeurs dans la liste
               insertion(listeBMU,i,j);
            }
        }
    }
    bmuF *bmu_g = listeBMU->premier;
    if(listeBMU->nombre_bmu>1){
    random = rand()%listeBMU->nombre_bmu;
    for(int next=1;next<random;next++)
    {
        bmu_g=bmu_g->suivant;
    }
    
    }
return bmu_g;
}

/* Fonction qui permets l'apprentissage*/
void apprentissage(int nb_iteration,double alpha,int nb_ligne,int nb_ligneR,int nb_colonneR){
    bmuF *Bmu;
    int ite = ((int)(N_config.nb_ite1/3));
    init_tab(nb_ligne);
        for(int i=0;i<nb_iteration;i++){
        shuffle(nb_ligne);
            for(int k=0;k<nb_ligne;k++){
            reseau.capteur=tab_v[index_tab[k]].arr;
            Bmu= detecte_BMU(nb_ligneR,nb_colonneR);
            strcpy(reseau.map[Bmu->ligne][Bmu->colonne].id,tab_v[index_tab[k]].etiquette);
            voisinage(Bmu,alpha);

            }
             if(reseau.taille_voisinage>1){
                if(N_config.nb_iteration == N_config.nb_ite1){
                    reseau.taille_voisinage-=1;
                    N_config.nb_ite1 = ite +((int)(500/3));
                }
                }

        }

}


int main(){

int Nb_lihgne = N_config.nombre_l=6;
int Nb_colonne = N_config.nombre_c=10;
double alpha1 = N_config.alphaM=0.7;
double alpha2 =  N_config.alphaM=0.07;
int Nb_iteration = N_config.nb_iteration = 2000;
int phase1 = N_config.nb_ite1 = N_config.nb_iteration/4;
int phase2 = Nb_iteration - phase1;
int longueurF = N_config.taille_f = 150;
int longueurV = N_config.taille_v = 4;

alloue_tab(longueurF);
parser_fi(longueurF,longueurV);
vect_moyen(longueurF,longueurV);
Minv(longueurV);
Maxv(longueurV);
config_NeuronMap();
affiche_map();
printf("\n");
reseau.taille_voisinage = 3;
printf(" ----------------phase1----------------\n");
apprentissage(phase1,alpha1,longueurF,Nb_lihgne,Nb_colonne);
affiche_map();
reseau.taille_voisinage = 1;
printf(" ----------------phase2----------------\n");
apprentissage(phase2,alpha2,longueurF,Nb_lihgne,Nb_colonne);
affiche_map();
free(tab_v);
free(max);
free(min);



return 0;
}