#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>



void Calcul_Caree(mpz_t value, mpz_t var  ){  mpz_mul(value,var,var); }

void Calcul_Modulo(mpz_t value , mpz_t n){ mpz_mod(value,value,n); }

void Caclcul_Multi(mpz_t value1 , mpz_t value2 ){ mpz_mul(value1, value1, value2); }

// on lutilise pas mais on a implementer
// x a lexposant e 
void Caclcul_a_lexposant(mpz_t x , mpz_t e )
{
	int taille =  mpz_sizeinbase (e,2);
    mpz_t res , tmp ;
    mpz_init(res);
    mpz_init(tmp);

    mpz_set_ui(tmp,1);
     
    for (int i = taille ; i >= 0 ; --i)
    {
				Calcul_Caree(res,tmp);
				mpz_set(tmp,res);
		    if (mpz_tstbit (e,i)==1)  Caclcul_Multi(tmp,x);      			
	}

    //le resulatat dans a 
	mpz_set(x,tmp);
}

// a exposant e mod n le resulatat est stocker dans a  
void Expo_Rapide(mpz_t a ,mpz_t e, mpz_t n){
     // on recupere la taille en binaire de lexposant pour pouvoir repeter 
	// la procedure 
    int taille =  mpz_sizeinbase (e,2);
    mpz_t res , tmp ;
    mpz_init(res);
    mpz_init(tmp);

    //on est obligé d'initialiser a 1 pour faire les calcul sinon on aura 0 
    mpz_set_ui(tmp,1);
     
    for (int i = taille ; i >= 0 ; --i)
    {
									
				// cas standard on calcule le carre du bits 			
				Calcul_Caree(res,tmp);
				mpz_set(tmp,res);

				//si le bit a lindice i  de lexpasant est 1 on le multiplie en plus du carre  
		     if (mpz_tstbit (e,i)==1)  Caclcul_Multi(tmp,a);
		        
		        //calcul du modulo
		         Calcul_Modulo(tmp,n);
										
	}

    //le resulatat dans a 
	mpz_set(a,tmp);
	
	
}

//res = pgcd (a,b)
void Calcul_PGCD(mpz_t res , mpz_t a , mpz_t b )
{

	mpz_t  reste ;
    mpz_init(reste);
      while (mpz_cmp_ui (b,0))
     {
    	mpz_mod(reste,a,b);
    	mpz_set(a,b);
    	mpz_set(b,reste);

     }
    
    mpz_set(res,a);
    mpz_clear(reste);

}

// test de fermat
int Test_Fermat(mpz_t nb_a_tester , int k  )
{
	mpz_t  a  ,  n_1 , tmp, UN ;
	mpz_init(tmp);
	mpz_init(UN);
    mpz_init(a);
    mpz_init(n_1);
    gmp_randstate_t rand;
   
   	mpz_set_ui(UN,1);
    gmp_randinit_default (rand);

   	// cas de base 2 et 1 
    if (!mpz_cmp_ui(nb_a_tester,2)){  printf("le nombre est premier!\n\n"); return 1;}
    if (!mpz_cmp_ui (nb_a_tester,1)){ printf("le nombre 1 est  non premier !\n\n"); return 0;}
	
  for (int i = 0 ; i< k ; i++ )
	{
		
		 mpz_urandomm (a , rand ,nb_a_tester);
		 // pour avoir de   1<a<n-1 
		 while (!mpz_cmp_ui (a,1) || !mpz_cmp_ui (a,0) ) {mpz_urandomm (a , rand ,nb_a_tester);}
        //calcul du n-1 
		mpz_sub (n_1,nb_a_tester,UN); 
		  // pour eviter de perdre la valeur initial de a
		  mpz_set(tmp,a); 
	    Expo_Rapide (a,n_1,nb_a_tester);
	    if (mpz_cmp_ui (a,1)) {printf("le nombre est composé \n\n"); return  0 ;}
	    mpz_set(a,tmp);
	}
	
	gmp_randclear(rand);
	mpz_clear(a);
	mpz_clear(n_1);
 	mpz_clear(tmp);
 	mpz_clear(UN);
	printf("le nombre est premier \n\n");
	return 1;
}





void Decomposition(mpz_t n_1 ,mpz_t t , mpz_t s )
{
    int cpt = 0;
  while (mpz_even_p(n_1)) 
  {
  		  mpz_div_ui (n_1, n_1, 2);
        cpt += 1;
            
  }
    mpz_set_ui(s,cpt);
    mpz_set(t,n_1);
   // gmp_printf("%Zd  %Zd\n" ,m ,k);
}

int Test_Miller_Rabin(mpz_t nb_a_tester , int k )
{
		
		mpz_t  n_1 ,  s , t ,a ,y ;
		mpz_init(y);
    mpz_init(a);
    mpz_init(s);
    mpz_init(t);
    mpz_init(n_1);
    	
      gmp_randstate_t rand;
   		gmp_randinit_default (rand);
   		mpz_sub_ui(n_1,nb_a_tester,1);
   		Decomposition(n_1,t,s);

   		if (!mpz_cmp_ui(nb_a_tester,2)){  printf("le nombre est premier!\n\n"); return 1;}
      if (!mpz_cmp_ui (nb_a_tester,1)){ printf("le nombre 1 est  non premier !\n\n"); return 0;}
    	 
       mpz_urandomm (a , rand ,nb_a_tester);
    	
      for (int i =1 ; i< k ; ++i )
    	{
          while ( !mpz_cmp_ui (a,0)) mpz_urandomm (a , rand ,nb_a_tester);
		      
          mpz_sub_ui(n_1,nb_a_tester,1);
          Expo_Rapide(a,t,nb_a_tester); 

		      mpz_set(y,a);

		 if ( mpz_cmp_ui (y,1) && mpz_cmp(y,n_1) )
		 {
		       int  etat = 0 ; 
            for ( int j = 0; j < (mpz_get_ui(s)) ; ++j)
            {
                 
                  //gmp_printf("%Zd  %Zd\n" ,m ,k);
				          Calcul_Caree( y , y);
				          Calcul_Modulo(y,nb_a_tester);
		            	 if ( !mpz_cmp_ui (y,1)){printf("le nombre est composé  \n" ); return 0;}
                   if ( !mpz_cmp_ui (y,mpz_get_ui(nb_a_tester)-1)) { etat = 1; ;break;}
            }
           if (!etat) { printf("nombre est composé \n");  return 0;}
      }
	   }
    	 printf ("le nombre est premier! \n"); return 1; 
}




int main(int argc , char ** argv ){


 if (argc < 3 ){ printf ("Erreur dans le nombre d'argument ! voir Makefile \n");}

 int k = atoi(argv[2]);

 mpz_t Nombre_a_tester  ; 

 mpz_init_set_str(Nombre_a_tester, argv[1], 10);



 printf("\nl'entier testé n = %s \npar Fermat avec K = %d \n\n",argv[1] , k);
 printf ("Resultat :   ");
 int a = Test_Fermat(Nombre_a_tester,k);

 printf("\nl'entier testé n = %s \npar Miller_rabin avec K = %d \n\n",argv[1] , k);
 printf ("Resultat :   ");
  a = Test_Miller_Rabin(Nombre_a_tester,k);



 
mpz_clear(Nombre_a_tester);








 a++ ; // eviter le warning 
 

 


  return 0; 
}