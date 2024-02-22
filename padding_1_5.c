#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include <time.h>
#include <unistd.h>
#include <string.h>

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
    		printf("\nQuel est le nom du fichier dans lequel vous désirez stocker le chiffré?\n\n");									// #
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
	srand(time(NULL));																	// #  On attribue à n la valeur de la clef publique et on calcule la taille de n en base 256
	int taille_n = taille_256(n)-11;													// #
	unsigned int taille_fichier = 0;													// ##


	while (fgetc(clair)!= EOF)															// ##
    {																					// #
    	taille_fichier += 1;															// #  On récupère la taille du fichier à chiffrer ou signer
    }																					// #
    rewind(clair);																		// ##


    int taille_padding = 8;																// #  On initialise les variables pour le padding
    mpz_set_ui(m,16);																	// #

    if (taille_fichier-1<taille_n)														// ##
    {																					// #
    	taille_padding = 8 + taille_n-(taille_fichier-1);								// #
    }																					// #
    																					// #
    for (int i = 0; i <taille_padding ;i++)												// #  Génération du padding pour le premier sous-message
    {																					// #
    	mpz_mul_ui(m,m,256);															// #
    	int A = (int)(rand() / (double)RAND_MAX * (255 - 17))+16;						// #
    	mpz_add_ui(m,m,A);																// #
    }																					// ##


    mpz_mul_ui(m,m,256);																// ##
    mpz_add_ui(m,m,0);																	// #  Initialisation de variable pour le chiffrement ou la signature
    unsigned int compteur = 0;															// ##


    for (int i = 0; i < taille_fichier; i++)											// #  Boucle pour lire tout le fichier
    {
    	mpz_mul_ui(m,m,256);															// ##
    	mpz_add_ui(m,m,fgetc(clair));													// #  Lecture d'un octet du fichier à chiffrer ou signer
    	compteur = compteur+1;															// ##

    	if ((compteur == taille_n) || (i==taille_fichier-1))							// #  Condition permettant de rentrer lorsque l'on a un sous-message m complet ou que l'on se trouve à la fin du fichier
    	{																				
    		exp_mod(m,m,e,n);															// ##
    		//gmp_fprintf(cypher,"%Zd ",m);												// #  Calcul et écriture d'un sous-message du chiffré ou de la signature dans le fichier destination
    		mpz_out_raw(cypher,m);														// ##
    		
    		mpz_set_ui(m,16);															// ##
    		if (taille_fichier-1-i<taille_n)											// #
    		{																			// #
    			taille_padding = 8 + taille_n-(taille_fichier-1-i);						// #
    		}																			// #
    		for (int j = 0; j <taille_padding;j++)										// #  Génération du padding pour le prochain sous-message
		    {																			// #
		    	mpz_mul_ui(m,m,256);													// #
		    	int A = (int)(rand() / (double)RAND_MAX * (255 - 17))+16;				// #
		    	mpz_add_ui(m,m,A);														// #
		    }																			// ##
		    
		    mpz_mul_ui(m,m,256);														// ##
    		mpz_add_ui(m,m,0);															// #
    		compteur = 0;																// #  Initialisation des variables pour le prochain tour
    	}																				// #
    }																					// ##

    if(signature == 1)																	// ##
    {																					// #
    	fclose(privee);																	// #
    }																					// #  Fermeture des fichiers
    mpz_clears(n, m, e, NULL);															// #
    fclose(clair);																		// #
	fclose(cypher);																		// #
	fclose(publique);																	// #
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
			printf("\nQuel est le nom du fichier dans lequel vous désirez stocker le clair? \n\n");											// #
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
		clair = fopen(nom_fichier_dechiffrer,"wb+");																						// #
	}																																		// #
	else																																	// #
	{																																		// #
		clair = fopen("signature_a_verifier","wb+");																						// #
	}																																		// ##

	while (fgetc(cypher)!= EOF)														// ##
    {																				// #
    	taille_fichier ++;															// #  Récupération de la taille du fichier à déchiffrer ou vérifier
    };																				// #
    rewind(cypher);																	// ##

	taille_n = taille_256(n)-2;														// ##
    ui_expo_ui(puissance,256,taille_n);												// # Initialisation des variables pour la récupération des octets
	mpz_set_ui(chiffre,0);															// ##

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

		mpz_fdiv_q(div,chiffre,puissance);											// ##
		do   																		// #
		{																			// #
			mpz_mul(puissance,puissance,div);										// #
			mpz_sub(chiffre,chiffre,puissance);										// #  Suppression du padding
			mpz_divexact(puissance,puissance,div);									// #
			mpz_divexact_ui(puissance,puissance,256);								// #
			mpz_fdiv_q(div,chiffre,puissance);										// #
		}while(mpz_cmp_ui(div,0)!= 0);												// #
		mpz_divexact_ui(puissance,puissance,256);									// ##

		while(mpz_cmp_ui(puissance,0) != 0)											// ##
		{																			// #
			mpz_fdiv_q(div,chiffre,puissance);										// #
			lettre = mpz_get_ui(div) ;												// #
			fputc(lettre,clair);													// #
																					// #
			if(mpz_cmp_ui(div,0) != 0)												// #  Récupération et écriture des octets de donnés dans le fichier destination
			{																		// #
				mpz_mul(puissance,puissance,div);									// #
				mpz_sub(chiffre,chiffre,puissance);									// #
				mpz_divexact(puissance,puissance,div);								// #
			}																		// #
			mpz_divexact_ui(puissance,puissance,256);								// #
		}																			// ##

		ui_expo_ui(puissance,256,taille_n);											// #  Réinitialisation de la variable puissance
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
				
				printf("\nQuelle est la longeur de la clé publique souhaitée?\n\n");
				scanf(" %d", &nombre_bit);


				generation_cle(nombre_bit, generateur);


				printf("\nQue souhaitez-vous faire maintenant?\n\n1 : Générer de nouvelles clé RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");
				goto marqueur;
				break;


			case 2:		// #  Chiffrement d'un fichier


				encrypt(0);


				printf("\nQue souhaitez-vous faire maintenant?\n\n1 : Générer de nouvelles clé RSA\n2 : Chiffrer un fichier\n3 : Déchiffrer un fichier\n4 : Signer un fichier\n5 : Vérifier une signature\n6 : Arret du programme\n\n");
				goto marqueur;
				break;



			case 3:		// #  Déchiffrement d'un fichier
				
				printf("\nComment souhaitez-vous déchiffrer?\n\n1 : Mode classique\n2 : Mode CRT\n\n");																									// ##
				scanf(" %d", &choix2);																																									// #
				while((choix2 != 1) & (choix2 != 2))																																					// #
				{																																														// #  Choix du mode de déchiffrement
					printf("\nLe chiffre que vous avez inqiqué ne correspond à aucune commande, comment souhaitez-vous déchiffrer?\n\n1 : Mode classique\n2 : Mode CRT\n3 : Arret du programme\n\n");	// #
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