#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define K 200        // Le nombre de premiers utilisés pour le crible K<=1000
#define securite 10  // t le paramètre entier de sécurité pour Miller-Rabin


// Cette fonction renvoie le reste modulaire d'un unsigned int
// Entrée : deux entiers a et n
// Sortie : un entier a modulo n
unsigned int mod (unsigned int a, unsigned int n)
{												  
	unsigned int q;
	q = (int) a / n; // Calcul du quotient q
	return a - q*n;
};


// Cette fonction calcule le reste modulaire d'un mpz
// Entrée : trois mpz reponse, nombre et module
// Sortie : vide mais reponse = nombre [module]
void modulo (mpz_t reponse, mpz_t nombre, mpz_t module)
{
	mpz_fdiv_r(reponse, nombre, module); // reponse = nombre [module]
};


// Même fonction que précdemment mais renvoie le résultat dans un int
// Entrée : un mpz nombre et un entier module
// Sortie : un entier reponse = nombre [module]
int modulo_ui (mpz_t nombre, int module)
{
	int reponse = 0;					  // ##
	mpz_t z_reponse, z_module;            // #  Initialisation des variables
	mpz_inits(z_reponse, z_module, NULL); // ##

	mpz_set_ui(z_reponse,reponse);        // ##
	mpz_set_ui(z_module,module);		  // #  Calcul du résultat via la fonction précédente modulo()
	modulo(z_reponse,nombre,z_module);	  // ##
	
	reponse = mpz_get_ui(z_reponse);      // Convertion du résultat de mpz à entier

	mpz_clears(z_reponse,z_module, NULL);
	return reponse;
};


// Cette fonction effectue une exponentiation d'entier et l'affetcte à un mpz
// Entrée : un mpz resultat et deux entiers a et b
// Sortie : vide mais resultat = a^b
void ui_expo_ui(mpz_t resultat, unsigned int a, unsigned int b)
{
	mpz_t m;                               // ##
	mpz_t d;							   // #  Initialisation des variables
	mpz_inits(m,d,NULL);				   // ##

	mpz_set_ui(m,a);					   // ##
	mpz_set_ui(d,b);					   // #  Convertion en mpz
	mpz_set_ui(resultat,1);				   // ##

	while(mpz_cmp_ui(d,0) != 0)            // ##
	{									   // #
		if(modulo_ui(d,2) == 1)			   // #
		{								   // #
			mpz_mul(resultat,resultat,m);  // #  Calcul de a^b par l'exponentation rapide
			mpz_sub_ui(d,d,1);			   // #  
		}								   // #
		mpz_mul(m,m,m);					   // #
		mpz_divexact_ui(d,d,2);			   // #
	}									   // ##

	mpz_clears(m,d,NULL);
};


// Cette fonction calcule une exponentiation modulaire et l'affecte à un mpz
// Entrée : quatre mpz resultat,m,d et n
// Sortie : void mais resultat = m^d [n]
void exp_mod(mpz_t resultat, mpz_t m, mpz_t d, mpz_t n) 
{
	mpz_t a;								// ##
	mpz_t b;								// #
	mpz_inits(a,b,NULL);					// #
											// #  Initialisation des variables
	mpz_set(a,m);							// #
	mpz_set(b,d);							// #
	mpz_set_ui(resultat,1);					// ##

	while(mpz_cmp_ui(b,0) != 0)				// ##
	{										// #
		if(modulo_ui(b,2) == 1)				// #
		{									// #
			mpz_mul(resultat,resultat,a);	// #
			modulo(resultat,resultat,n);	// #  Calcul de a^b par l'exponentiation rapide avec calcule de restes modulo n
			mpz_sub_ui(b,b,1);				// #
		}									// #
		mpz_mul(a,a,a);						// #
		modulo(a,a,n);						// #
		mpz_divexact_ui(b,b,2);				// #
	}										// ##

	mpz_clears(a,b,NULL);
};


// Cette fonction est la même que la précédente mais l'exposant est un entier
// Entrée : trois mpz resultat, m et n et un entier d
// Sortie : vide mais resultat = m^d [n]
void exp_mod_ui(mpz_t resultat, mpz_t m, unsigned int d, mpz_t n)
{
	mpz_t z_d;
	mpz_init_set_ui(z_d,d);		// Convertion de d d'entier à mpz

	exp_mod(resultat,m,z_d,n);	// Calcul du résultat via la fonction précédente exp_mod()

	mpz_clear(z_d);
};


// Cette fonction calcule PGCD(a,b) et les coefficients de Bezout associés avec euclide étendu
// Entrée : cinq mpz a,b,x,y et pgcd
// Sortie : vide mais pgcd = PGCD(a,b) et pgcd = ax + by
void euclid_gcd(mpz_t a,mpz_t b, mpz_t x, mpz_t y, mpz_t pgcd)
{
	mpz_t u,v,q,r,m,n;				// ##
	mpz_inits(u,v,q,r,m,n, NULL);	// #
	mpz_set_ui(u,1);				// #  Initialisation des variables
	mpz_set_ui(y,1);				// #
	mpz_set_ui(x,0);				// #
	mpz_set_ui(v,0);				// ##

	while (mpz_cmp_ui(b,0) !=0)		// ##
	{								// #
		mpz_fdiv_qr(q,r,a,b);		// #
		mpz_mul(m,u,q);				// #
		mpz_mul(n,v,q);				// #
		mpz_sub(m,x,m);				// #
		mpz_sub(n,y,n);				// #  
		mpz_set(a,b);				// #  Calcul du PGCD par euclide étendu
		mpz_set(b,r);				// #
		mpz_set(x,u);				// #
		mpz_set(y,v);				// #
		mpz_set(u,m);				// #
		mpz_set(v,n);				// #
	}								// #
	mpz_set(pgcd,a);				// ##

	mpz_clears(u,v,q,r,m,n,NULL);		
};


// Cette fonction calcule l'inverse modulaire si il existe
// Entrée : trois mpz inv_mod,b et a
// Sortie : vide mais inv_mod = b^(-1) [a]
void modular_inv(mpz_t inv_mod, mpz_t b, mpz_t a)
{
	mpz_t save_a,x,y,pgcd,save_b;				// ##
	mpz_inits(save_a,x,y,pgcd,save_b,NULL);		// #  Initialisation des variables et sauvegarde de a et b
	mpz_set(save_a,a);							// #
	mpz_set(save_b,b);							// ##
	
	euclid_gcd(a,b,x,y,pgcd);					// Calcul du PGCD et des coefficients de Bezout
	
	if (mpz_cmp_ui(pgcd,1) != 0)				// ##
	{											// #
		mpz_set_ui(inv_mod,0);					// #  Vérification de PGCD(a,b) = 1
		mpz_clears(x,y,pgcd,NULL);				// #
		return;									// #
	}											// ##
	
	modulo(inv_mod , x , save_a);				// ##
	mpz_set(a,save_a);							// #  Réduction et affectattion de b^(-1) [a] et remise de a et b à leur valeur initiale
	mpz_set(b,save_b);							// ##
	
	mpz_clears(x,y,pgcd,save_a,save_b,NULL);
};


// Cette fonction vérifie si un nombre et premier par Miller-Rabin
// Entrée : un mpz nombre et un générateur aléatoire generateur
// Sortie : 0 si nombre est composé et 1 s'il est probablement premier
int Miller_Rabin(mpz_t nombre, gmp_randstate_t generateur)
{
	mpz_t r, a, y, j;
	mpz_inits(r, a, y, j, NULL);											// ##
	mpz_set(r, nombre);														// #  Initialisation des variables
	mpz_sub_ui(r,r,1);														// #
	int s = 0;																// ##
	
	while(modulo_ui(r,2) == 0)												// ##
	{																		// #
		s++;																// #  Calcul de r et s tel que nombre-1 = r x 2^s
		mpz_divexact_ui(r,r,2);												// #
	}																		// ##

	for(unsigned int i =1; i< securite; i++)								// Boucle avec le paramètre de sécurité
	{
		mpz_sub_ui(nombre,nombre,3);										// ##
		mpz_urandomm(a, generateur, nombre);								// #  Tirage aléatoire de a entre 2 et n-2
		mpz_add_ui(nombre,nombre,3);										// #
		mpz_add_ui(a, a, 2);												// ##

		exp_mod(y,a,r,nombre);												// Calcul de y = a^r [nombre]
		
		mpz_sub_ui(nombre,nombre,1); 										
		if((mpz_cmp_ui(y,1) != 0) && (mpz_cmp(y,nombre) != 0))				
		{
			mpz_set_ui(j,1);
			while((mpz_cmp_ui(j,(s-1)) <= 0) && (mpz_cmp(y,nombre) != 0))	// ##
			{																// #
				mpz_add_ui(nombre,nombre,1);								// #
				exp_mod_ui(y,y,2,nombre);									// #
				mpz_sub_ui(nombre,nombre,1); 								// #
				if(mpz_cmp_ui(y,1) == 0)									// #
				{															// #
					mpz_add_ui(nombre,nombre,1); 							// #
					mpz_clears(r, a, y, j, NULL);							// #
					return 0;												// #
				}															// #  Test de Miller-Rabin
				mpz_add_ui(j,j,1);											// #
			}																// #									
			if(mpz_cmp(y,nombre) != 0)										// #
			{																// #
				mpz_add_ui(nombre,nombre,1); 								// #
				mpz_clears(r, a, y, j, NULL);								// #
				return 0;													// #
			}																// #
		}																	// #
		mpz_add_ui(nombre,nombre,1);										// ##
	}

	mpz_clears(r, a, y, j, NULL);
	return 1;
};


// Cette fonction génère aléatoirement un nombre premier par la méthode du crible optimisé
// Entrée : un mpz nb_premier, un entier b et un générateur alétoire state
// Sortie : vide mais nb_premier est un nombre premier de taille b bits
void optimized_crible_generation(mpz_t nb_premier, unsigned int b, gmp_randstate_t state)
{
	mpz_t sub;												// ##																			
	mpz_t s_1;												// #  Initialisation des variables
	mpz_t s_2;												// #
	mpz_inits(sub,s_1,s_2,NULL);							// ##

	unsigned int res[K];									// ##
	unsigned int r;											// #  Déclaration d'un tableau de résidus et d'un tableau de premiers de taille K
	unsigned int premiers[K];								// ##
	
	FILE* fichier;											// #  Ouverture du fichier contenant les XX premiers nombres premiers
	fichier = fopen("primes.txt","r");						// #
	
	for(unsigned int i = 0; i< K; i++)						// ##
	{														// #
		fscanf(fichier,"%d",premiers+i);					// #  Lecture et stockage des K plus petits premiers
		fgetc(fichier);										// #
	}														// ##
	
	ui_expo_ui(s_1,2,b);									// ##
	ui_expo_ui(s_2,2,b-1);									// #  Calcul de sub = 2^b - 2^(b-1) -1
	mpz_sub(sub,s_1,s_2);									// #
	mpz_sub_ui(sub,sub,1);									// ##
	
	do 														// ##
	{														// #
		mpz_urandomm(nb_premier,state,sub);					// #  Génération d'un nombre aléatoire impaire s'écrivant sur b bits
		mpz_add(nb_premier,nb_premier,s_2);					// #
		r = modulo_ui(nb_premier,2);						// #
	}while(r == 0);											// ##
	
	for(int i=0;i<K;i++) 									// ## 
	{														// #
		r = modulo_ui(nb_premier,premiers[i]);				// #  Initialisation du tableau de résidus
		res[i] = r;											// #
	}														// ##

	etiquette:
		for(unsigned int j=0;j<K;j++)						// ##
		{													// #
			while(res[j] == 0)								// #
			{												// #
				for(unsigned int l=0;l<K;l++)				// #
				{											// #  Modification de nb_premier jusqu'à trouver un candidat
					res[l] = mod(res[l]+2, premiers[l]);	// #
				}											// #
				mpz_add_ui(nb_premier,nb_premier,2);		// #
			}												// #
		}													// ##

		if(Miller_Rabin(nb_premier,state) == 0)				// ##
		{													// #
			for(int m=0;m<K;m++)							// #
			{												// #
				res[m] = mod(res[m]+2, premiers[m]);		// #  Appel de la fonction Miller_Rabin() pour tester notre candidat, on recommence au cas échéant
			}												// #
			mpz_add_ui(nb_premier,nb_premier,2);			// #
			goto etiquette;									// #
		}													// ##
		goto etiquette2;

	etiquette2:
		fclose(fichier);
		mpz_clears(sub,s_1,s_2,NULL);
};


// Cette fonction génère l'ensemble des clefs publique et secrète nécessaire pour RSA
// Entrée : un entier nombre_bit et un génnérateur aléatoire generateur
// Sortie : vide mais création d'un fichier contenant une clef publique et un autre la clef secrète associée
void generation_cle(unsigned int nombre_bit, gmp_randstate_t generateur)
{

	mpz_t borne, phi, e, cle_publique, cle_prive, p, q, Ip;							// ##
	mpz_inits(borne,phi,e, cle_publique, cle_prive, p, q, Ip, NULL);				// #  Initialisation des variables
	mpz_set_ui(e,65537);															// #
	char choix;																		// ##

	char nom_fichier_cle_publique[100];																									// ##
	etiquette:																															// #
		printf("\nQuel est le nom du fichier dans lequel vous désirez stocker la clé publique?\n\n");									// #
		scanf(" %99s", nom_fichier_cle_publique);																							// #
		if(access( nom_fichier_cle_publique, F_OK ) == 0)																				// #
		{																																// #
			printf("\nAttention, le nom de fichier saisi existe déjà, êtes-vous sûr de vouloir l'effacer?[Y/N]\n\n");					// #
			scanf(" %c",&choix);																										// #
			while((choix != 'Y') & (choix != 'N'))																						// #
			{																															// #  Création et ouverture du fichier qui contiendra la clef publique
				printf("Le choix que vous avez fait n'a pas été compris, êtes-vous sûr de vouloir l'effacer (attendu Y ou N)?\n\n");	// #
				scanf(" %c",&choix);																									// #
			}																															// #
			if(choix == 'N')																											// #
			{																															// #
				goto etiquette;																											// #
			}																															// #
		}																																// #
	FILE* publique = fopen(nom_fichier_cle_publique,"wb+");																				// ##


	char nom_fichier_cle_secrete[100];																									// ##
	etiquette2:																															// #
    	printf("\nQuel est le nom du fichier dans lequel vous désirez stocker la clé privée?\n\n");										// #
   		scanf(" %99s", nom_fichier_cle_secrete);																							// #
   		if(access( nom_fichier_cle_secrete, F_OK ) == 0)																				// #
		{																																// #
			printf("\nAttention, le nom de fichier saisi existe déjà, êtes-vous sûr de vouloir l'effacer?[Y/N]\n\n");					// #
			scanf(" %c",&choix);																										// #
			while((choix != 'Y') & (choix != 'N'))																						// #
			{																															// #  Création et ouverture du fichier qui contiendra la clef secrète
				printf("Le choix que vous avez fait n'a pas été compris, êtes-vous sûr de vouloir l'effacer (attendu Y ou N)?\n\n");	// #
				scanf(" %c",&choix);																									// #
			}																															// #
			if(choix == 'N')																											// #
			{																															// #
				goto etiquette2;																										// #
			}																															// #
		}																																// #
	FILE* secret = fopen(nom_fichier_cle_secrete,"wb+");																				// ##


	ui_expo_ui(borne,2,nombre_bit-1);												// ##
	do 																				// #
	{																				// #
		if( mod(nombre_bit,2) == 0)													// #
		{																			// #
			optimized_crible_generation(p, nombre_bit/2, generateur);				// #
			optimized_crible_generation(q, nombre_bit/2, generateur);				// #
		}																			// #
		else																		// # Appel de la fonction optimized_crible_generation() pour génèrer p et q et calcul de la clef publique avec vériffication de sa taille
		{																			// #
			optimized_crible_generation(p, (nombre_bit-1)/2, generateur);			// #
			optimized_crible_generation(q, (nombre_bit+1)/2, generateur);			// #
		}																			// #
		mpz_mul(cle_publique, p, q);												// #
	}																				// #
	while(mpz_cmp(cle_publique,borne) >= 0);										// ##

	mpz_sub_ui(p,p,1);																// ## 
	mpz_sub_ui(q,q,1);																// #
	mpz_mul(phi,p,q);																// #  Calcul de phi
	mpz_add_ui(p,p,1);																// #
	mpz_add_ui(q,q,1);																// ##

	modular_inv(cle_prive, e, phi);													// #  Calcul de la clef secrète
	modular_inv(Ip,p,q);															// #

	mpz_out_raw(publique,cle_publique);												// #  Ecriture des clefs dans leur fichier respectifs
	mpz_out_raw(secret,cle_prive);													// #
	mpz_out_raw(secret,p);
	mpz_out_raw(secret,q);
	mpz_out_raw(secret,Ip);
	
	mpz_clears(borne, phi, e, cle_publique, cle_prive, p, q, Ip, NULL);
	fclose(publique);
	fclose(secret);
};


// Cette fonction donne la taille d'un nombre en base 256
// Entrée : un mpz n
// Sortie : un entier compteur donnant en combien d'octet s'écrit n
unsigned int taille_256(mpz_t n)
{
	mpz_t a;						// ##
	mpz_init_set_ui(a,256);			// #  Initialisation des variables
	unsigned int compteur = 1;		// ##

	while(mpz_cmp(n,a) > 0)			// ##
	{								// #
		mpz_mul_ui(a,a,256);		// #  Boucle comparant n avec les puissances de 256 jusqu'à ce que n < 256^a
		compteur = compteur+1;		// #
	}								// ##

	mpz_clear(a);
	return(compteur);
};


// Cette fonction hashe un fichier à partir de son nom
// Entrée : une chaine de caractere contenant le nom du fichier à hasher
// Sortie : vide mais crée un fichier de 256 bits contenant le hasher du fichier d'entrée sous le nom HASHER
void SHA256(char* nom_du_fichier_a_hasher)
{
	mpz_t h;										// ##
	mpz_init(h);									// #
	char resultat[64];								// #	Initialisation des variables
	char commande[100] = "openssl sha256 ";			// #
	strcat(commande,nom_du_fichier_a_hasher);		// ##

	FILE* hash = popen(commande,"r");				// ##
	fscanf(hash, "%64s",resultat);					// #	utilisation et récupération de la commande openssl sha256 "nom_fichier"
	fscanf(hash, "%64s",resultat);					// ##

	
	mpz_set_str(h,resultat,16);						// ##
	FILE* fichier_hash = fopen("HASHER","wb+");		// #	stockage du résultat dans le fichier HASHER
	mpz_out_raw(fichier_hash,h);					// ##

	mpz_clear(h);									// ##
	pclose(hash);									// #	Fermeture des fichiers et libération de l'espace
	fclose(fichier_hash);							// ##
};


// Cette fonction permet de hasher une chaîne de caractère et la stocke dans un fichier
// Entrée : une chaîne de caractère à hasher
// Sortie : vide mais ajoute le hasher de la chaine de caractère à la fin du fichier MGF1
void sha256sum(char* txt)
{
	mpz_t h;										// ##
	mpz_init(h);									// #
	char resultat[64];								// #	
	char commande1[100] = "echo -n ";				// #  Initialisation des variables et construction de la commande à envoyer au terminale
	char commande2[100] =" | sha256sum";			// #
	strcat(commande1,txt);							// #
	strcat(commande1,commande2);					// ##

	FILE* hash = popen(commande1,"r");				// #  Utilisation et récupération de la commande echo -n "txt" | sha256sum
	fscanf(hash, "%64s",resultat);					// #

	mpz_set_str(h,resultat,16);						// ##
	FILE* T = fopen("MGF1","a");					// #  Stockage du résultat dans le fichier HASHER
	gmp_fprintf(T, "%Zx",h );						// ##

	mpz_clear(h);									// ##
	pclose(hash);									// #  Fermeture des fichiers et libération de l'espace
	fclose(T);										// ##
}


// Cette fonction convertit un caractère en hexadécimal en un entier en décimal. 
// Entrée : un caractère c écrit en hexadécimal
// Sortie : un entier correspondant à la valeur en décimal de c
int hex_to_int(char c)
{
	if(c >= '0' && c <= '9')						// ##	
	{												// #
		return(c - '0');							// #
	}												// #  Conversion de c en décimal 
	if(c >= 'a' && c <= 'f')						// #
	{												// #
		return(c - 'a' + 10);						// #
	}												// ##
}


// Cette fonction convertit un entier en décimal en un caractère en hexadécimal. 
// Entrée : un entier a compris entre 0 et 15
// Sortie : un caractère correspondant à la valeur en hexadécimal de a
char int_to_hex(int a)
{
	if(a < 10)										// ##
	{												// #
		return(a + '0');							// #
	}												// #  Conversion de a en hexadécimal
	else											// #
	{												// #
		return(a + 'a' - 10);						// #
	}												// ##
}


// Cette fonction correspond à la fonction utiliser dans MGF1
// Entrée : un mpz_t counter qui correspond à une seed et un pointeur sur un caractère pour récupérer le résultat
// Sortie : vide mais le résultat est stocké dans la chaine resultat
void I2OSP(mpz_t counter,char *resultat)
{
	mpz_t div;											// ##
	mpz_t puissance;									// #
	mpz_t save;											// #
	mpz_inits(div,puissance,save,NULL);					// #  Initialisation des variables
	ui_expo_ui(puissance,16,7);							// #
	mpz_set(save,counter);								// #
	int i = 0;											// #
	int j;												// ##

	while(mpz_cmp_ui(puissance,0) != 0)					// ##
	{													// #
		mpz_fdiv_q(div,save,puissance);					// #
		j = mpz_get_ui(div);							// #
		*(resultat+i) = int_to_hex(j);					// #
														// #
		if(mpz_cmp_ui(div,0) != 0)						// #
		{												// #  Algorithme I2OSP
			mpz_mul(puissance,puissance,div);			// #
			mpz_sub(save,save,puissance);				// #
			mpz_divexact(puissance,puissance,div);		// #
		}												// #
		mpz_divexact_ui(puissance,puissance,16);		// #
		i = i+1;										// #
	}													// ##

	mpz_clears(div,puissance,save,NULL);
}


// Cette fonction correspond à la fonction MGF1 utilisée pour OAEP
// Entrée : un mpz_t correspondant à la taille de la sortie souhaitée
// Sortie : vide mais la sortie de l'algorithme est stockée dans le fichier "MGF1"
void MGF1(mpz_t seed, mpz_t l)
{
	mpz_t a;											// ##
	mpz_init(a);										// #  Initialisation des variables
	ui_expo_ui(a,2,37);									// ##

	if(mpz_cmp(l,a) > 0)								// ##
	{													// #  Vérification de la taille l
		printf("\nMasque trop long\n");					// #
	}													// ##
	else
	{
		mpz_t counter_max;								// ##
		mpz_t i;										// #
		mpz_inits(counter_max,i,NULL);					// #
		mpz_fdiv_q_ui(counter_max,l,32);				// #
		mpz_add_ui(counter_max,counter_max,1);			// #
														// #
		char chaine[24] = { 0 };						// #
		char * pointeur;								// #
		mpz_get_str(chaine,16,seed);					// #
		pointeur = &chaine[16];							// #
														// #  Algorithme MGF1
		mpz_set_ui(i,0);								// #
		while(mpz_cmp_ui(counter_max,0) > 0)			// #
		{												// #
			I2OSP(i,pointeur);							// #
			sha256sum(chaine);							// #
			mpz_sub_ui(counter_max,counter_max,1);		// #
			mpz_add_ui(i,i,1);							// #
		}												// ##
		mpz_clears(counter_max,NULL);
	}
	mpz_clear(a);
}


// Cette fonction permet d'effectuer le XOR entre 2 valeurs hexadécimales et renvoie le résultat en hexadécimal.
// Entrée : deux caractères en hexadécimal a et b a xorer entre eux
// Sortie : un caractère en hexadécimal c = XOR(a,b)
char XOR(char a,char b)
{
	char c;													// ##
	int int_a = 0;											// #
	int int_b = 0;											// #
	int int_c = 0;											// #  Déclaration et initialisation des variables
	int div_a;												// #
	int div_b;												// #
	int puissance = 8;										// ##

	int_a = hex_to_int(a);									// #  Conversion de a et b en décimal
	int_b = hex_to_int(b);									// #

	while(puissance != 0)									// ##
	{														// #
		div_a = int_a / puissance;							// #
		div_b = int_b / puissance;							// #
		int_c = int_c + (div_a^div_b) * puissance;			// #
															// #
		if(div_a != 0)										// #
		{													// #
			puissance = puissance * div_a;					// #
			int_a = int_a - puissance;						// #
			puissance = puissance / div_a;					// #  Calcul du XOR entre a et b
		}													// #
															// #
		if(div_b != 0)										// #
		{													// #
			puissance = puissance * div_b;					// #
			int_b = int_b - puissance;						// #
			puissance = puissance / div_b;					// #
		}													// #
															// #
		puissance = puissance / 2;							// #
	}														// ##

	c = int_to_hex(int_c);									// #  Conversion de c en hexadécimal

	return(c);
}


// Cette fonction applique la padding OAEP à un sous-message et l'écrit dans un fichier
// Entrée : un entier length_n correspondant à la taille du sous-message sans le padding, un entier dernier nous indiquant si nous sommes au dernier bloc et un flux correspondant au fichier à chiffrer
// Sortie : vide mais un sous-message est ajouté à la suite du fichier OAEP
void OAEP(int length_n,int dernier, FILE* fichier)
{
	srand(time(NULL));													// #  Initialisation du random pour générer l'aléa

	mpz_t l;															// ##
	mpz_t n;															// #
	mpz_t m;															// #
	mpz_t puissance;													// #
	mpz_t div;															// #
	mpz_init_set_ui(l,length_n);										// #
	mpz_init_set_ui(n,0);												// #
	mpz_init_set_ui(m,0);												// #
	mpz_init_set_ui(puissance,0);										// #  Déclaration des variables
	mpz_init_set_ui(div,0);												// #
	char c1;															// #
	char c2;															// #
	char c3;															// #
	char chaine[16] = { 0 };											// #
	int A;																// #
	int h1;																// #
	int h2;																// ##

	if (dernier==0)														// ##
	{																	// #
		for (int i = 0; i <8 ;i++)										// #
	    {																// #
	    	mpz_mul_ui(m,m,256);										// #  Test si nous sommes au dernier bloc à chiffrer :
	    	A = (int)(rand() / (double)RAND_MAX * (255 - 17))+16;		// #
	    	mpz_add_ui(m,m,A);											// #  	- si non on génère 8 octets aléatoirement dans m
	    }																// #
	}																	// #	- si oui on affecte la taille du dernier bloc dans m
	else																// #
	{																	// #
		mpz_set_ui(m,dernier);											// #
	}																	// ##

	MGF1(m,l);															// #  Passage de m dans MGF1 avec l pour taille de sortie

	FILE* MGF = fopen("MGF1","r");										// #  Ouverture des fichiers OAEP et MGF1
	FILE* OAE = fopen("OAEP","a");										// #

	for(int a=0;a<length_n;a++)											// #  Boucle pour créer X (sachant que le bloc créé à la fin est de la forme X||Y)
	{
		if(!feof(fichier))												// ##
		{																// #
			A = fgetc(fichier);											// #
		}																// #
		else															// #  Lecture d'un octet A du fichier à chiffrer et initialisation de h1 et h2 tel que A = h1*16 + h2
		{																// #
			A = 0;														// #
		}																// #
		h1 = A/16;														// #
		h2 = A - (h1*16);												// ##

		c1 = fgetc(MGF);												// #  Lecture d'une valeur hexa du fichier MGF1

		c2 = int_to_hex(h1);											// #  XOR entre le message et MGF1 qui donne les octets de X
		c3 = XOR(c1,c2);												// #

		if(a<8)															// ##
		{																// #
			mpz_mul_ui(n,n,16);											// #  Stockage des 8 premiers octets de X (1)
			A = hex_to_int(c3);											// #
    		mpz_add_ui(n,n,A);											// #
		}																// ##

		fprintf(OAE,"%c",c3);											// #  Ecriture dans le fichier OAEP

		c1 = fgetc(MGF);												// #  Lecture d'une valeur hexa du fichier MGF1
		
		c2 = int_to_hex(h2);											// #  XOR entre le message et MGF1 qui donne les octets de X
		c3 = XOR(c1,c2);												// #

		if (a<8)														// ##
		{																// #
			mpz_mul_ui(n,n,16);											// #  Stockage des 8 premiers octets de X (2)
    		A = hex_to_int(c3);											// #
    		mpz_add_ui(n,n,A);											// #
		}																// ##

		fprintf(OAE,"%c",c3);											// #  Ecriture dans le fichier OAEP
	}

	fclose(MGF);														// #  Fermeture et suppression du fichier MGF1
	remove("MGF1");														// #

	mpz_set_ui(l,8);													// ##
	MGF1(n,l);															// #  Passage de X dans MGF1 avec 8 pour taille de sortie et ouverture de MGF1
	MGF = fopen("MGF1","r");											// ##

	if(dernier == 0)													// ##
	{																	// #
		mpz_get_str(chaine,16,m);										// #
	}																	// #
	else																// #
	{																	// #
		for (int j=15;j>=0;j--)											// #
		{																// #
			ui_expo_ui(puissance,16,j);									// #  Initialisation de chaine qui va contenir le padding
			mpz_fdiv_q(div,m,puissance);								// #
																		// #
			A = mpz_get_ui(div);										// #
			chaine[15-j] = int_to_hex(A);								// #
																		// #
			mpz_mul(puissance,puissance,div);							// #
			mpz_sub(m,m,puissance);										// #
		}																// #
	}																	// ##

 	for(int b=0;b<16;b++)												// ##
	{																	// #
		c1 = chaine[b];													// #
		c2 = fgetc(MGF);												// #  Calcul et écriture de Y = XOR(MGF1(X),padding) dans le fichier OAEP à la suite de X
		c3 = XOR(c1,c2);												// #
		fprintf(OAE,"%c",c3);											// #
	}																	// ##

	fclose(MGF);														// #  Fermeture et suppression du fichier MGF1
	remove("MGF1");														// #
	
	fclose(OAE);														// #  Fermeture du fichier OAEP qui contient le bloc X||Y à chiffrer
	
	mpz_clears(l,n,m,puissance,div,NULL);
}


// Cette fonction supprime le padding d'un bloc et écrit le résultat à la suite du fichier clair
// Entrée : un entier length_n correspondant à la taille d'un bloc, un entier dernier nous indiquant si nous sommes au dernier bloc, un mpz correspondant à un bloc et un flux correspondant au fichier destination
// Sortie : vide mais le sous-message est écrit à la fin du fichier clair
void inv_OAEP(int length_n, int dernier, mpz_t chiffre,FILE* clair)
{
	mpz_t puissance;									// ##
	mpz_t z_X;											// #
	mpz_t div;											// #
	mpz_t temp;											// #
	mpz_t r;											// #
	mpz_t l;											// #
	mpz_inits(puissance,z_X,div,temp,r,l,NULL);			// #
	mpz_set_ui(r,0);									// #
	mpz_set_ui(div,0);									// #
	mpz_set_ui(temp,0);									// #
	mpz_set_ui(l,8);									// #  Déclaration des variables
	char MGF1_X;										// #
	char Y;												// #
	char c1;											// #
	char c2;											// #
	char c3;											// #
	unsigned int lettre = 0;							// #
	int A;												// #
	int B;												// #
	int int_r;											// #
	int rtemp;											// ##
	
	ui_expo_ui(puissance,256,(length_n-8));				// #  Récupération des 8 premiers octets de X
	mpz_fdiv_q(z_X,chiffre,puissance);					// #
	
	MGF1(z_X,l);										// #  Passage des 8 premiers octets de X dans MGF1 avec 8 pour taille de sortie et ouverture de MGF1
	FILE* MGF = fopen("MGF1","r");						// #

	ui_expo_ui(puissance,256,8);						// ##
	mpz_fdiv_q(div,chiffre,puissance);					// #  Récupération de X 
	mpz_set(z_X,div);									// ##

	mpz_mul(puissance,puissance,div);					// #  Récupération de Y
	mpz_sub(temp,chiffre,puissance);					// #

	for(int i=15;i>=0;i--)								// #  Boucle pour calculer XOR(MGF1(X),Y)
	{
		MGF1_X = fgetc(MGF);							// #  Lecture d'une valeur hexadécimale de MGF1(X)

		ui_expo_ui(puissance,16,i);						// ##
		mpz_fdiv_q(div,temp,puissance);					// #  Récupération d'une valeur hexadécimale de Y
		int_r = mpz_get_ui(div);						// #
		Y = int_to_hex(int_r);							// ##
		
		c1 = XOR(MGF1_X,Y);								// #  XOR entre MGF1(X) et Y

		A = hex_to_int(c1);								// ##
		mpz_mul_ui(r,r,16);								// #  stockage du résultat dans r
		mpz_add_ui(r,r,A);								// ##

		mpz_mul(puissance,puissance,div);				// #  Décrémentation de Y
		mpz_sub(temp,temp,puissance);					// #
	}

	fclose(MGF);										// #  Fermeture et suppression du fichier MGF1
	remove("MGF1");										// #

	mpz_set_ui(l,(length_n - 8));						// ##
	MGF1(r,l);											// #  Passage du padding (r) dans MGF1 avec la taille d'un bloc moins la taille du padding pour taille de sortie et ouverture de MGF1
	MGF = fopen("MGF1","r");							// ##

	if(dernier != 0)									// ##
	{													// #  Si on est au dernier bloc on récupère la taille du dernier sous-message
		rtemp = mpz_get_ui(r);							// #
	}													// ##

	ui_expo_ui(puissance,16,((2*(length_n-8))-1));		// #  Initialisation de puissance

	for(int j=0;j<(length_n - 8);j++)					// #  Boucle pour le dernier XOR entre X et MGF1(r)
	{
		mpz_fdiv_q(div,z_X,puissance);					// ##
		A = mpz_get_ui(div);							// #  Récupération d'une valeur hexadécimale de X
		c2 = int_to_hex(A);								// ##

		c1 = fgetc(MGF);								// #  Récupération d'une valeur hexadécimale de MGF1(r)

		c3 = XOR(c1,c2);								// #  XOR entre X et MGF1(r)

		if(mpz_cmp_ui(div,0)!=0)						// ##
		{												// #
			mpz_mul(puissance,puissance,div);			// #
			mpz_sub(z_X,z_X,puissance);					// #  Décrémentation de puissance et de X
			mpz_divexact(puissance,puissance,div);		// #
		}												// #
		mpz_divexact_ui(puissance,puissance,16);		// ##
		
		B = hex_to_int(c3);								// #  Stockage du résultat dans lettre
		lettre = B * 16;								// #
		
		mpz_fdiv_q(div,z_X,puissance);					// ##
		A = mpz_get_ui(div);							// #  Récupération d'une valeur hexadécimale de X
		c2 = int_to_hex(A);								// ##
		
		c1 = fgetc(MGF);								// #  Récupération d'une valeur hexadécimale de MGF1(r)

		c3 = XOR(c1,c2);								// #  XOR entre X et MGF1(r)

		if(mpz_cmp_ui(div,0)!=0)						// ##
		{												// #
			mpz_mul(puissance,puissance,div);			// #
			mpz_sub(z_X,z_X,puissance);					// #  Décrémentation de puissance et de X
			mpz_divexact(puissance,puissance,div);		// #
		}												// #
		mpz_divexact_ui(puissance,puissance,16);		// ##

		B = hex_to_int(c3);								// #  Stockage du résultat dans lettre
		lettre = lettre + B;							// #

		if (dernier==0)									// ##
		{												// #
			fputc(lettre,clair);						// #
		}												// #
		else											// #
		{												// #  Ecriture d'un octet du résultat dans le fichier clair
			if (j < rtemp)								// #
			{											// #
				fputc(lettre,clair);					// #
			}											// #
		}												// #
		lettre = 0;										// ##
	}

	mpz_clears(puissance,z_X,div,temp,r,l,NULL);		// ##
	fclose(MGF);										// #  Fermeture et suppression du fichier MGF1
	remove("MGF1");										// ##
}


// Cette fonction sert soit à chiffrer un fichier, soit à signer un fichier 
// Entrée : un entier signature, si signature vaut 0 on chiffre, si signature vaut 1 on signe
// Sortie : vide mais on crée soit un fichier contenant un chiffré ou un fichier contenant une signature
void encrypt(unsigned int signature) 
{
	char choix;

	char nom_fichier_cle_publique[100];													// ##
																						// #
	etiquette:																			// #
   		printf("\nQuel est le nom du fichier contenant la clé publique?\n\n");			// #
    	scanf(" %99s", nom_fichier_cle_publique);										// #
		if(access( nom_fichier_cle_publique, F_OK ) != 0)								// #  Ouverture du fichier contenant la clé publique
		{																				// #
			printf("\nAttention, le nom de fichier saisi n'existe pas!");				// #
			goto etiquette;																// #
		}																				// #
    FILE* publique = fopen(nom_fichier_cle_publique,"rb+");								// ##


	char nom_fichier_a_chiffrer[100];													// ##
																						// #
	etiquette2:																			// #
		if(signature == 0)																// #
		{																				// #
   			printf("\nQuel est le nom du fichier que vous désirez chiffrer?\n\n");		// #
   		}																				// #
   		else																			// #
   		{																				// #  Ouverture du fichier contenant le fichier que l'on va signer ou chiffrer
   			printf("\nQuel est le nom du fichier que vous désirez signer?\n\n");		// #
   		}																				// #
    	scanf(" %99s", nom_fichier_a_chiffrer);											// #
    	if(access( nom_fichier_a_chiffrer, F_OK ) != 0)									// #
		{																				// #
			printf("\nAttention, le nom de fichier saisi n'existe pas!");				// #
			goto etiquette2;															// #
		}																				// #
	FILE* clair;																		// #
	if(signature == 0)																	// #
	{																					// #
   		clair = fopen(nom_fichier_a_chiffrer,"rb+");									// #
   	}																					// #
   	else																				// #
   	{																					// #
   		SHA256(nom_fichier_a_chiffrer);													// #
   		clair = fopen("HASHER","rb+");													// #
   	}																					// #
																						// ##


    char nom_fichier_chiffrer[100];																										// ##
    																																	// #
    etiquette3:																															// #
    	if (signature == 0)																												// #
    	{																																// #
    		printf("\nQuel est le nom du fichier dans lequel vous désirez stocker le chiffré ?\n\n");									// #
    	}																																// #
    	else																															// #
    	{																																// #
    		printf("\nQuel est le nom du fichier dans lequel vous désirez stocker la signature?\n\n");									// #
    	}																																// #
   		scanf(" %99s", nom_fichier_chiffrer);																							// #
   		if(access( nom_fichier_chiffrer, F_OK ) == 0)																					// #
		{																																// #  Création ou édition et ouverture du fichier de destination
			printf("\nAttention, le nom de fichier saisi existe déjà, êtes-vous sûr de vouloir l'effacer?[Y/N]\n\n");					// #
			scanf(" %c",&choix);																										// #
			while((choix != 'Y') & (choix != 'N'))																						// #
			{																															// #
				printf("Le choix que vous avez fait n'a pas été compris, êtes-vous sûr de vouloir l'effacer (attendu Y ou N)?\n\n");	// #
				scanf(" %c",&choix);																									// #
			}																															// #
																																		// #
			if(choix == 'N')																											// #
			{																															// #
				goto etiquette3;																										// #
			}																															// #
		}																																// #
	FILE* cypher = fopen(nom_fichier_chiffrer,"wb+");																					// ##



	mpz_t n, m, e;																		// ##
	mpz_inits(n, m, e, NULL);															// #  Déclaration des variables pour chiffrer ou signer
	FILE* privee;																		// ##
	
	if(signature == 0)																	// ##
	{																					// #
		mpz_set_ui(e,65537);															// #
	}																					// #
	else																				// #
	{																					// #
		char nom_fichier_cle_privee[100];												// #  
																						// #  Choix entre chiffrement et signature :
		etiquette4:																		// #  
			printf("\nQuel est le nom du fichier contenant la clé privée?\n\n");		// #  	- chiffrement : on initialise e à 65537
			scanf(" %99s", nom_fichier_cle_privee);										// #  
			if(access( nom_fichier_cle_privee, F_OK ) != 0)								// #  	- signature : on ouvre le fichier contenant la clef privée et on initialise e à la valeur de la clef privée
			{																			// #
				printf("\nAttention, le nom de fichier saisi n'existe pas!");			// #
				goto etiquette4;														// #
			}																			// #
		privee = fopen(nom_fichier_cle_privee,"rb+");									// #
		//gmp_fscanf(privee, " %Zd", e);												// #
		mpz_inp_raw(e,privee);															// #
	}																					// ##

	//gmp_fscanf(publique, " %Zd", n);													// ##
	mpz_inp_raw(n,publique);															// #
																						// #  On attribue à n la valeur de la clef publique et on calcule la taille de n en base 256
	int taille_n = taille_256(n)-1;														// #
	unsigned int taille_fichier = 0;													// ##
	int length_n;

	fseek(clair,0,SEEK_END);															// ##
    taille_fichier = ftell(clair);														// #  On récupère la taille du fichier à chiffrer ou signer
    rewind(clair);																		// ##

    int dernier = 0;																	// ##
    int taille_padding = 8;																// #
    length_n = taille_n - taille_padding;												// #  On initialise les variables pour le padding
    int i = 0;																			// #
    unsigned int k = 0;																			// ##
    
    while(i < taille_fichier)															// #  Boucle pour lire tout le fichier
    {
    	if (taille_fichier-i<taille_n)													// ##
   		{																				// #  Modification pour le dernier bloc
    		dernier = taille_fichier - i;												// #
    	}																				// ##
    	
    	OAEP(length_n,dernier,clair);													// #  Padding
    	i = i + length_n;
    }		

    FILE* OAE = fopen("OAEP","r");														// #  Ouverture du fichier contenant les blocs à chiffrer
    
    taille_fichier = 0;																	// ##
    i = 0;																				// #  Initialisation des variables pour chiffrer/signer
    char c;																				// ##

	fseek(OAE,0,SEEK_END);																// ##
    taille_fichier = ftell(OAE);														// #  On récupère la taille du fichier à chiffrer ou signer
    rewind(OAE);																		// ##

    while(i < taille_fichier)															// #  Boucle pour lire tout le fichier
    {	
    	mpz_set_ui(m,0);																// ##
    	for (int j = 0; j < 2*taille_n; j++)											// #
    	{																				// #
    		c = fgetc(OAE);																// #
    		mpz_mul_ui(m,m,16);															// #
			if(c >= '0' && c <= '9')													// #
			{																			// #  Lecture d'un bloc à chiffrer
				mpz_add_ui(m,m,c - '0');												// #
			}																			// #
			if(c >= 'a' && c <= 'f')													// #
			{																			// #
				mpz_add_ui(m,m, c - 'a' + 10);											// #
			}																			// #
    	}																				// ##
    	
	    exp_mod(m,m,e,n);																// #  Chiffrement d'un bloc
	    mpz_out_raw(cypher,m);															// #

	    i += 2*taille_n;																			
	}

    if(signature == 1)																	// ##
    {																					// #
    	fclose(privee);																	// #
    }																					// #
    mpz_clears(n, m, e, NULL);															// #
    fclose(OAE);																		// #
    fclose(clair);																		// #  Fermeture et suppression des fichiers
	fclose(cypher);																		// #
	fclose(publique);																	// #
	remove("OAEP");																		// #
	if(signature == 1)																	// #
	{																					// #
		remove("HASHER");																// #
	}																					// ##
};


// Cette fonction sert à déchiffrer un fichier en mode standard ou en mode crt ou à déchiffrer une signature
// Entrée : deux entiers crt et signature, si signature vaut 0 on déchiffre, si signature vaut 1 on déchiffre une signature. Si crt vaut 0 on dechiffre en mode standard si crt vaut 1 on déchiffre en mode crt.
// Sortie : vide mais on crée un fichier contenant le clair ou la signature déchiffrée
void decrypt(unsigned int crt, unsigned int signature)
{														
	char choix;																				// ##
	mpz_t n, d, p, q, Ip, chiffre, puissance, div, dp, dq, mp, mq ;							// #
    mpz_inits(n, d, p, q, Ip, chiffre, puissance, div,dp, dq, mp, mq,NULL);					// #
    unsigned int lettre;																	// #  Initialisation des variables
	unsigned int taille_n;																	// #
	unsigned int compteur = 0;																// #
	unsigned int taille_fichier = 0;														// #
	FILE* privee;																			// ##

	char nom_fichier_a_dechiffrer[100];														// ##
	etiquette:																				// #
		if(signature == 0)																	// #
		{																					// #
			printf("\nQuel est le nom du fichier que vous souhaitez déchiffrer? \n\n");		// #
		}																					// #
		else																				// #
		{																					// #
			printf("\nQuel est le nom du fichier contenant la signature? \n\n");			// #  Ouverture du fichier à déchiffrer ou de la signature à déchiffrer
		}																					// #
		scanf(" %99s", nom_fichier_a_dechiffrer);											// #
		if(access( nom_fichier_a_dechiffrer, F_OK ) != 0)									// #
		{																					// #
			printf("\nAttention, le nom de fichier saisi n'existe pas!");					// #
			goto etiquette;																	// #
		}																					// #
	FILE* cypher = fopen(nom_fichier_a_dechiffrer,"rb+");									// ##



    char nom_fichier_cle_publique[100];														// ##
    etiquette2:																				// #
		printf("\nQuel est le nom du fichier contenant la clé publique?\n\n");				// #
		scanf(" %99s", nom_fichier_cle_publique);											// #
		if(access( nom_fichier_cle_publique, F_OK ) != 0)									// #  Ouverture du fichier contenant la clef publique
		{																					// #
			printf("\nAttention, le nom de fichier saisi n'existe pas!");					// #
			goto etiquette2;																// #
		}																					// #
	FILE* publique = fopen(nom_fichier_cle_publique,"rb+");									// ##
	

    mpz_inp_raw(n,publique);																// #  On initialise n à la valeur de la clef publique


	if(signature == 0)																		// ##
	{																						// #
		char nom_fichier_cle_privee[100];													// #
	    etiquette3:																			// #
			printf("\nQuel est le nom du fichier contenant la clé privée?\n\n");			// #
			scanf(" %99s", nom_fichier_cle_privee);											// #
			if(access( nom_fichier_cle_privee, F_OK ) != 0)									// #
			{																				// #
				printf("\nAttention, le nom de fichier saisi n'existe pas!");				// #  Choix entre déchiffrement d'un fichier chiffré ou d'une signature :
				goto etiquette3;															// #
			}																				// #  	- fichier chiffré : on récupère les clefs privées stockées dans un fichier et on les initialise (d vaut la valeur de la clef privée)
		privee = fopen(nom_fichier_cle_privee,"rb+");										// #
																							// #  	- signature chiffrée : on initialise d à la valeur de l'exposant publique 65537
	    mpz_inp_raw(d,privee);																// #
	    mpz_inp_raw(p,privee);																// #
	    mpz_inp_raw(q,privee);																// #
	    mpz_inp_raw(Ip,privee);																// #
	}																						// #
	else																					// #
	{																						// #
		mpz_set_ui(d,65537);																// #
	}																						// ##


	FILE* clair;																															// ##
	if(signature == 0)																														// #
	{																																		// #
		char nom_fichier_dechiffrer[100];																									// #
		etiquette4:																															// #
			printf("\nQuel est le nom du fichier dans lequel vous désirez stocké le clair? \n\n");											// #
			scanf(" %99s", nom_fichier_dechiffrer);																							// #
			if(access( nom_fichier_dechiffrer, F_OK ) == 0)																					// #
			{																																// #
				printf("\nAttention, le nom de fichier saisi existe déjà, êtes-vous sûr de vouloir l'effacer?[Y/N]\n\n");					// #
				scanf(" %c",&choix);																										// #
				while((choix != 'Y') & (choix != 'N'))																						// #
				{																															// #
					printf("Le choix que vous avez fait n'a pas été compris, êtes-vous sûr de vouloir l'effacer (attendu Y ou N)?\n\n");	// #  Création et ouverture du fichier destination où l'on stocke le clair ou la signature à vérifier après déchiffrement
					scanf(" %c",&choix);																									// #
				}																															// #
																																			// #
				if(choix == 'N')																											// #
				{																															// #
					goto etiquette4;																										// #
				}																															// #
			}																																// #
		clair = fopen(nom_fichier_dechiffrer,"wb");																						// #
	}																																		// #
	else																																	// #
	{																																		// #
		clair = fopen("signature_a_verifier","wb+");																						// #
	}																																		// ##

	fseek(cypher,0,SEEK_END);														// ##
    taille_fichier = ftell(cypher);													// #  On récupère la taille du fichier à déchiffrer ou de la signature à déchiffrer
    rewind(cypher);																	// ##

	taille_n = taille_256(n)-1;														// ##
    ui_expo_ui(puissance,256,taille_n);												// # Initialisation des variables pour la récupération des octets
	mpz_set_ui(chiffre,0);
	int dernier = 0;															// ##

	if(crt == 1)																	// ##
	{																				// #
		mpz_sub_ui(p,p,1);															// #
		modulo(dp,d,p);																// #
		mpz_add_ui(p,p,1);															// #  Initialisation des variables pour le déchiffrement en mode crt
																					// #
		mpz_sub_ui(q,q,1);															// #
		modulo(dq,d,q);																// #
		mpz_add_ui(q,q,1);															// #
	}																				// ##

	do 																				// #  Boucle pour lire tout le fichier
	{
		compteur = compteur + mpz_inp_raw(chiffre,cypher);							// #  Mise à jour de la valeur de compteur et lecture du mpz_t chiffré
		if(crt == 0)																// ##
		{																			// #  Déchiffrement en mode standard ou déchiffrement de la signature
			exp_mod(chiffre,chiffre,d,n);											// #
		}																			// ##

		if(crt == 1)																// ##
		{																			// #
			exp_mod(mp,chiffre,dp,p);												// #
			exp_mod(mq,chiffre,dq,q);												// #
			mpz_set(chiffre,mq);													// #
			mpz_sub(chiffre,chiffre,mp);											// #  Déchiffrement en mode crt
			mpz_mul(chiffre,chiffre,Ip);											// #
			modulo(chiffre,chiffre,q);												// #
			mpz_mul(chiffre,chiffre,p);												// #
			mpz_add(chiffre,chiffre,mp);											// #
		}																			// ##

		if (compteur==taille_fichier)
		{
			dernier=1;
		}
		inv_OAEP(taille_n,dernier,chiffre,clair);
	}while(compteur < taille_fichier);												

	fclose(cypher);																	// ##
	fclose(clair);																	// #
	fclose(publique);																// #
	if(signature == 0)																// #  Fermeture des fichiers
	{																				// #
		fclose(privee);																// #
	}																				// #
	mpz_clears(n, d, p, q, Ip, chiffre, puissance, div, dp, dq, mp, mq, NULL);		// ##
};


// Cette fonction sert à vérifier une signature
// Entrée : vide
// Sortie : vide mais affichage de la validité de la signature
void verification_signature()
{
	decrypt(0,1);																	// #  Déchiffrement de la signature

	int caractere1, caractere2;																						// ##
	char nom_fichier_a_verifier[100];																				// #
	etiquette:																										// #
		printf("\nQuel est le nom du fichier original dont vous voulez vérifier la signature associée? \n\n");		// #
		scanf(" %99s", nom_fichier_a_verifier);																		// #
		if(access( nom_fichier_a_verifier, F_OK ) != 0)																// #  Ouverture du fichier original et du fichier contenant la signature déchiffrée
		{																											// #
			printf("\nAttention, le nom de fichier saisi n'existe pas!");											// #
			goto etiquette;																							// #
		}																											// #
	SHA256(nom_fichier_a_verifier);																					// #
	FILE* cypher = fopen("HASHER","rb+");																			// #
	FILE* signature = fopen("signature_a_verifier","rb+");															// ##

	caractere1 = fgetc(cypher);														// ##
	caractere2 = fgetc(signature);													// #
	while((caractere1 != EOF) & (caractere2 != EOF) & (caractere1 == caractere2))	// #
	{																				// #  Comparaison octet par octet des deux fichiers
		caractere1 = fgetc(cypher);													// #
		caractere2 = fgetc(signature);												// #
	}																				// ##

	if(caractere1 == caractere2)													// ##
	{																				// #
		printf("\nLa signature est valide!\n\n");									// #
	}																				// #  Affichage du résultat
	else																			// #
	{																				// #
		printf("\nLa signature est invalide.\n\n");									// #
	}																				// ##
	
	fclose(cypher);																	// ##
	fclose(signature);																// #  Fermeture des fichiers et suppression du fichier intermédiaire
	remove("signature_a_verifier");													// #
	remove("HASHER");																// ##											
};


// Corps du programme
int main()
{
	gmp_randstate_t generateur;					// ##
	gmp_randinit_default(generateur);			// #  Initialisation des variable et du générateur aléatoire
	gmp_randseed_ui(generateur, time(NULL));	// #
	unsigned int choix1, choix2, nombre_bit;	// ##

	printf("\nBienvenue dans le programme de chiffrement RSA, que souhaitez-vous faire?\n\n1 : Générer de nouvelles clés RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");	// #  Affichage des options du programme
	

	marqueur:
		scanf(" %d", &choix1);																																																																// ##
		while((choix1 != 1) & (choix1 != 2) & (choix1 != 3) & (choix1 != 4) & (choix1 != 5) & (choix1 != 6))																																												// #
		{																																																																					// #  Boucle pour avoir une réponse valide
			printf("\nLe chiffre que vous avez indiqué ne correspond à aucune commande, que souhaiez-vous faire?\n\n1 : Générer de nouvelles clé RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");	// #
			scanf(" %d", &choix1);																																																															// #
		}																																																																					// ##
		

		switch(choix1)	// #  Switch séparant toutes les options du programmes
		{
			case 1:		// #  Générations de nouvelles clefs RSA
				
				printf("\nQuelle est la longeur de la clé publique souhaitée (en bit) ?\n\n");
				scanf(" %d", &nombre_bit);


				generation_cle(nombre_bit, generateur);


				printf("\nQue souhaiez-vous faire maintenant?\n\n1 : Générer de nouvelles clé RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");
				goto marqueur;
				break;


			case 2:		// #  Chiffrement d'un fichier


				encrypt(0);


				printf("\nQue souhaiez-vous faire maintenant?\n\n1 : Générer de nouvelles clé RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");
				goto marqueur;
				break;



			case 3:		// #  Déchiffrement d'un fichier
				
				printf("\nComment souhaitez-vous déchiffrer?\n\n1 : Mode classique\n2 : Mode CRT\n\n");																									// ##
				scanf(" %d", &choix2);																																									// #
				while((choix2 != 1) & (choix2 != 2))																																					// #
				{																																														// #  Choix du mode de déchiffrement
					printf("\nLe chiffre que vous avez indiqué ne correspond à aucune commande, comment souhaitez-vous déchiffrer?\n\n1 : Mode classique\n2 : Mode CRT\n3 : Arret du programme\n\n");	// #
					scanf(" %d", &choix2);																																								// #
				}																																														// ##


				if(choix2 == 3)
				{
					printf("\nArret du programme.\n");
					break;
				}



				decrypt(choix2-1,0);



				printf("\nQue souhaitez-vous faire maintenant?\n\n1 : Générer de nouvelles clé RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");
				goto marqueur;
				break;

			case 4:		// #  Signature d'un fichier


				encrypt(1);

				printf("\nQue souhaitez-vous faire maintenant?\n\n1 : Générer de nouvelles clé RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");
				goto marqueur;
				break;

			case 5:		// #  Vérification d'une signature


				verification_signature();

				printf("\nQue souhaitez-vous faire maintenant?\n\n1 : Générer de nouvelles clé RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");
				goto marqueur;
				break;

			case 6:		// #  Arrêt du programme

				printf("\nArret du programme.\n");
				break;
		}

	gmp_randclear(generateur);
	return 0;
};